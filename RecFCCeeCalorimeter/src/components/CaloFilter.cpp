/*
 * Copyright (c) 2014-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "Gaudi/Property.h"

// edm4hep
#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Constants.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/TimeSeriesCollection.h"

#include "podio/UserDataCollection.h"

#include "k4FWCore/Transformer.h"

#include <ostream>
#include <string>

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include <TMatrixD.h>
#include <TRandom3.h>

using DigitsColl = edm4hep::TimeSeriesCollection;

using FilterColl = edm4hep::TimeSeriesCollection;
using MatchedSampleIdxColl = podio::UserDataCollection<int>;
using MatchedSampleEnergyColl = podio::UserDataCollection<float>;

struct CaloFilterFunc final
    : k4FWCore::MultiTransformer<std::tuple<FilterColl, MatchedSampleIdxColl, MatchedSampleEnergyColl>(
          const DigitsColl&)> {
  CaloFilterFunc(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc, {KeyValues("InputCollection", {"TimeSeriesCollection"})},

                         {KeyValues("OutputCollectionFilteredPulse", {"FilteredDigitsCollection"}),
                          KeyValues("OutputCollectionMatchedSampleIdx", {"MatchedSampleIdx"}),
                          KeyValues("OutputCollectionMatchedSampleEnergy", {"MatchedSampleEnergy"})}) {}

  StatusCode initialize() override {
    // Get matched filter template
    if (m_filterName.value() == "MatchedDirac") {
      info() << "Using the matched filter: Matched Dirac" << endmsg;
      FilterTemplate = {0.0, 1.0, 0.0};
    } else if (m_filterName.value() == "Matched_Gaussian") {
      info() << "Using the matched filter: Matched Gaussian" << endmsg;

      // Number of samples you want the known pulse (i.e. filter) to consider
      // Note that an odd number of samples > 1 is the best for a Gaussian.
      // Otherwise, migration effects will happen

      // If even number of samples then up count > down count
      int DownCnt = (m_filterTemplateSize.value() - 1) / 2;
      int UpCnt = m_filterTemplateSize.value() - DownCnt - 1;

      // We will create a Gaussian shape w/ mean = m_mu and sigma = m_sigma. The
      // filter will be defined as e^(-0.5 * ((x
      // - m_mu)/m_sigma)^2) for x in [maximum - DownCnt, maximum + UpCnt].

      int samplingInterval = (m_pulseEndTime.value() - m_pulseInitTime.value()) / m_lenSample.value();
      std::vector<float> PulseShape(m_lenSample.value(), 0.0);
      for (int i = 0; i < m_lenSample.value(); ++i) {
        PulseShape[i] = Gaussian(i * samplingInterval, m_mu.value(), m_sigma.value());
      }

      // Get the filter template between MaxIdx - DownCnt and MaxIdx + UpCnt
      auto MaxIdx = std::distance(PulseShape.begin(), std::max_element(PulseShape.begin(), PulseShape.end()));

      int DownIdx = MaxIdx - DownCnt;
      int UpIdx = MaxIdx + UpCnt;
      DI = DownIdx;
      UI = UpIdx;
      // Validate indices
      if (DownIdx < 0 || UpIdx >= static_cast<int>(PulseShape.size())) {
        error() << "Index out of bounds: DownIdx=" << DownIdx << ", UpIdx=" << UpIdx
                << ", PulseShape.size()=" << PulseShape.size() << endmsg;
        return StatusCode::FAILURE;
      }
      // Resize FilterTemplate to fit the range [DownIdx, UpIdx]
      FilterTemplate.resize(UpIdx - DownIdx + 1, 0.0);

      info() << "MaxIdx: " << MaxIdx << ", DownIdx: " << DownIdx << ", UpIdx: " << UpIdx << endmsg;
      // Get the filter template and calculate SumSquared for normalization

      for (int i = DownIdx; i <= UpIdx; ++i) {
        FilterTemplate[i - DownIdx] = PulseShape[i]; // Use relative indexing for FilterTemplate
        SumSquared += PulseShape[i] * PulseShape[i];
      }

      /*
      for (int i = DownIdx; i <= UpIdx; ++i) {
          FilterTemplate[i - DownIdx] = PulseShape[i]; // Use relative indexing
      for FilterTemplate SumSquared += PulseShape[i] * PulseShape[i];
      }
       *
       *
       * */
    } else {
      error() << "Unknown filter name: " << m_filterName.value() << endmsg;
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  // This is the function that will be called to transform the data
  // Note that the function has to be const, as well as all pointers to
  // collections we get from the input
  std::tuple<FilterColl, MatchedSampleIdxColl, MatchedSampleEnergyColl>
  operator()(const DigitsColl& DigitsPulse) const override {
    info() << "Digitized pulse collection size: " << DigitsPulse.size() << endmsg;

    DigitsColl FilteredDigitsCollection;
    MatchedSampleIdxColl MaxIdxCollection;
    MatchedSampleEnergyColl EnergyCollection;

    // Loop over DigitsPulse to extract the pulse amplitudes
    for (const auto& Digit : DigitsPulse) {
      const auto InputPulse = Digit.getAmplitude();

      auto FilteredDigit = FilteredDigitsCollection.create();
      FilteredDigit.setCellID(Digit.getCellID());

      FilteredDigit.setTime(0.0);                     // Placeholder for time info
      FilteredDigit.setInterval(Digit.getInterval()); // Set the interval for the digitized pulse in ns

      // Apply matched filter
      auto Out = applyMatchedFilter(InputPulse, FilterTemplate);

      for (unsigned int i = 0; i < Out.size(); i++) {

        FilteredDigit.addToAmplitude(Out[i]);
      }

      // Calculate the energy of the matched filter and the matched sample index
      auto MaxIdx = std::distance(Out.begin(), std::max_element(Out.begin(), Out.end()));
      auto MaxVal = *std::max_element(Out.begin(), Out.end());
      auto Energy = MaxVal / SumSquared;

      debug() << "Cell ID" << Digit.getCellID() << ", MaxIdx: " << MaxIdx << ", MaxVal: " << MaxVal
              << ", Energy: " << Energy << endmsg;

      // Store the matched sample index and energy
      MaxIdxCollection.push_back(MaxIdx);
      EnergyCollection.push_back(Energy);
    }

    return std::make_tuple(std::move(FilteredDigitsCollection), std::move(MaxIdxCollection),
                           std::move(EnergyCollection));
  }

  StatusCode finalize() override { return StatusCode::SUCCESS; }

private:
  /// Map to be used for the lookup of the pulse shapes
  Gaudi::Property<std::string> m_filterName{this, "filterName", "Matched_Gaussian", "Name of the filter to apply"};

  Gaudi::Property<float> m_mu{this, "mu", 50, "Mean of Gaussian pulse"};
  Gaudi::Property<float> m_sigma{this, "sigma", 20, "Sigma of Gaussian pulse"};

  // Number of samples in the filter template
  Gaudi::Property<int> m_filterTemplateSize{this, "filterTemplateSize", 5, "Size of the filter template"};

  // Initial time of the pulse
  Gaudi::Property<float> m_pulseInitTime{this, "pulseInitTime", 0.0, "Initial time of the pulse"};
  // End time of the pulse
  Gaudi::Property<float> m_pulseEndTime{this, "pulseEndTime", 750.0, "End time of the pulse"};
  // Number of samples in pulse
  Gaudi::Property<int> m_lenSample{this, "pulseSamplingLength", 30, "Number of samples in pulse"};
  // Noise Energy
  Gaudi::Property<float> m_noiseEnergy{this, "noiseEnergy", 0.001, "Noise Energy"};
  // Noise Width
  Gaudi::Property<float> m_noiseWidth{this, "noiseWidth", 1, "Noise Standard Deviation"};
  // Cell count
  Gaudi::Property<float> m_cellcount{this, "cellCount", 100, "Amount of Cells for simulated noise"};
  // DownIndex
  std::vector<float> FilterTemplate;
  int DI, UI;
  float SumSquared = 0;

  // Method to apply matched filter
  std::vector<float> applyMatchedFilter(const podio::RelationRange<float>& pulse,
                                        const std::vector<float>& Cfilter) const {
    // Reverse filter
    info() << "pulse.size() is (From applyMatchedFilter) " << pulse.size() << "\n";

    TMatrixD* m1 = new TMatrixD(pulse.size(), m_cellcount.value());

    TMatrixD* minv = new TMatrixD(pulse.size(), pulse.size());

    info() << "m1 is " << m1->GetNrows() << " x " << m1->GetNcols() << "\n";
    info() << "minv is " << minv->GetNrows() << " x " << minv->GetNcols() << "\n";
    AutoCorrMatFromCovariance(pulse, m1, minv, m_cellcount.value());

    std::vector<float> filter = FilterFunc(minv, pulse);
    std::vector<float> filter0 = Cfilter;

    std::reverse(filter.begin(), filter.end());
    std::reverse(filter0.begin(), filter0.end());

    float tsum[2];
    tsum[0] = 0.0f, tsum[1] = 0.0f;
    // info() << "h : { ";
    for (int i = 0; i < filter.size(); i++) {
      // info() << filter[i] << " ,";
      tsum[0] += filter[i];
    };
    tsum[0] /= ((float)filter.size());
    // info() << " }\n";

    info() << "O : { ";
    for (int i = 0; i < filter0.size(); i++) {
      // info() << filter0[i] << " ,";
      tsum[1] += filter0[i];
    };
    tsum[1] /= ((float)filter0.size());
    info() << "h = " << tsum[0] / tsum[1] << "O";

    // convolve the pulse with the filter
    return ConvolveSame(pulse, filter);
  }

  float Gaussian(float x, float mu = 0.0, float sigma = 1.0) const {
    return (1.0 / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((x - mu) / sigma, 2));
  }

  float GaussianNoNorm(float x, float mu = 0.0, float sigma = 1.0) const {
    return exp(-0.5 * pow((x - mu) / sigma, 2));
  }

  std::vector<float> ConvolveSame(const podio::RelationRange<float>& pulse, const std::vector<float>& filter) const {
    int n = pulse.size();
    int m = filter.size();
    int conv_size = n + m - 1;
    int start = (m - 1) / 2;
    int end = start + n;

    std::vector<float> result(conv_size, 0.0);
    for (size_t i = 0; i < pulse.size(); ++i) {
      for (size_t j = 0; j < filter.size(); ++j) {
        result[i + j] += pulse[i] * filter[j];
      }
    }
    std::vector<float> result_same(result.begin() + start, result.begin() + end);
    return result_same;
  }
  // By Edgar
  TRandom* r3 = new TRandom();

  void AutoCorrMatFromCovariance(const podio::RelationRange<float>& pulse, TMatrixD* m1, TMatrixD* minv,
                                 int Cellcount) const {
    int N = pulse.size();
    std::vector<float> h(N);
    float sum = 0;
    float mu = m_noiseEnergy.value();
    float sigma = m_noiseWidth.value();
    float mean[pulse.size()];
    float T[pulse.size()];
    for (int i = 0; i < m1->GetNrows(); i++) {
      for (int j = 0; j < m1->GetNcols(); j++) {
        N = r3->Gaus(mu, sigma) - mu;
        m1->SetElement(i, j, N);
        sum += N;
      };
      mean[i] = (sum / Cellcount);

      sum = 0;
    }
    // k is Cell index

    for (int i = 0; i < minv->GetNrows(); i++) {
      for (int j = i; j < minv->GetNcols(); j++) {
        for (int k = 0; k < Cellcount; k++) {
          sum += (((*m1)(i, k) - mean[i]) * ((*m1)(j, k) - mean[j]));
        }
        sum /= Cellcount;
        minv->SetElement(i, j, sum);
        minv->SetElement(j, i, sum);
        sum = 0;
      };
    }
    minv->Invert();
    minv->Print();
  }

  std::vector<float> FilterFunc(TMatrixD* minv, const podio::RelationRange<float>& DigitsPulse) const {
    std::vector<float> h(DigitsPulse.size());
    float squaresum = 0;
    info() << "DigisPulse.size() is \t " << DigitsPulse.size() << "\n";
    std::vector<float> hfilt(UI - DI + 1);
    for (int i = 0; i < DigitsPulse.size(); i++) {
      for (int j = 0; j < DigitsPulse.size(); j++) {
        h[i] += (*minv)(i, j) * DigitsPulse[j];
      }
      squaresum += h[i] * DigitsPulse[i];
    }
    /*
      for (int i = 0; i < DigitsPulse.size(); i++){
       h[i] = (i > DI && i <= UI ) *  h[i]/sqrt(squaresum);
      };*/

    for (int i = DI; i <= UI; i++) {
      hfilt[i - DI] = h[i] / sqrt(squaresum);
    };
    return hfilt;
  };
};

DECLARE_COMPONENT(CaloFilterFunc)
