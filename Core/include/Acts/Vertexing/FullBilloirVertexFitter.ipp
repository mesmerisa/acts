// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/VertexingError.hpp"


namespace {
/*
double getDeterminantBernie(const std::vector<std::vector<double>> vect) {
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    } 
    int dimension = vect.size();

    if(dimension == 0) {
        return 1;
    }

    if(dimension == 1) {
        return vect[0][0];
    }

    //Formula for 2x2-matrix
    if(dimension == 2) {
        return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
    }

    double result = 0;
    int sign = 1;
    for(int i = 0; i < dimension; i++) {

        //Submatrix
        std::vector<std::vector<double>> subVect(dimension - 1, std::vector<double> (dimension - 1));
        for(int m = 1; m < dimension; m++) {
            int z = 0;
            for(int n = 0; n < dimension; n++) {
                if(n != i) {
                    subVect[m-1][z] = vect[m][n];
                    z++;
                }
            }
        }

        //recursive call
        result = result + sign * vect[0][i] * getDeterminantBernie(subVect);
        sign = -sign;
    }

    return result;
}

std::vector<std::vector<double>> getTransposeBernie(const std::vector<std::vector<double>> matrix1) {

    //Transpose-matrix: height = width(matrix), width = height(matrix)
    std::vector<std::vector<double>> solution(matrix1[0].size(), std::vector<double> (matrix1.size()));

    //Filling solution-matrix
    for(size_t i = 0; i < matrix1.size(); i++) {
        for(size_t j = 0; j < matrix1[0].size(); j++) {
            solution[j][i] = matrix1[i][j];
        }
    }
    return solution;
}

std::vector<std::vector<double>> getCofactorBernie(const std::vector<std::vector<double>> vect) {
    if(vect.size() != vect[0].size()) {
        throw std::runtime_error("Matrix is not quadratic");
    } 

    std::vector<std::vector<double>> solution(vect.size(), std::vector<double> (vect.size()));
    std::vector<std::vector<double>> subVect(vect.size() - 1, std::vector<double> (vect.size() - 1));

    for(std::size_t i = 0; i < vect.size(); i++) {
        for(std::size_t j = 0; j < vect[0].size(); j++) {

            int p = 0;
            for(size_t x = 0; x < vect.size(); x++) {
                if(x == i) {
                    continue;
                }
                int q = 0;

                for(size_t y = 0; y < vect.size(); y++) {
                    if(y == j) {
                        continue;
                    }

                    subVect[p][q] = vect[x][y];
                    q++;
                }
                p++;
            }
            solution[i][j] = pow(-1, i + j) * getDeterminantBernie(subVect);
        }
    }
    return solution;
}


std::vector<std::vector<double>> getInverseBernie(const std::vector<std::vector<double>> vect) {
   std::cout << "Determinant: " << getDeterminantBernie(vect) << std::endl;
    if(getDeterminantBernie(vect) == 0) {
        std::cout << "Determinant is 0, det: " << getDeterminantBernie(vect) << std::endl;
        return std::vector<std::vector<double>>();
    } 

    double d = 1.0/getDeterminantBernie(vect);
    std::vector<std::vector<double>> solution(vect.size(), std::vector<double> (vect.size()));

    for(size_t i = 0; i < vect.size(); i++) {
        for(size_t j = 0; j < vect.size(); j++) {
            solution[i][j] = vect[i][j]; 
        }
    }

    solution = getTransposeBernie(getCofactorBernie(solution));

    for(size_t i = 0; i < vect.size(); i++) {
        for(size_t j = 0; j < vect.size(); j++) {
            solution[i][j] *= d;
        }
    }

    return solution;
}

void printMatrixBernie(const std::vector<std::vector<double>> vect) {
    for(std::size_t i = 0; i < vect.size(); i++) {
        for(std::size_t j = 0; j < vect[0].size(); j++) {
            std::cout << std::setw(8) << vect[i][j] << " ";
        }
        std::cout << "\n";
    }
}
*/

/// @struct BilloirTrack
///
/// @brief Struct to cache track-specific matrix operations in Billoir fitter
template <typename input_track_t>
struct BilloirTrack {
  using Jacobian = Acts::ActsMatrix<Acts::eBoundSize, 4>;

  BilloirTrack(const input_track_t* params, Acts::LinearizedTrack lTrack)
      : originalTrack(params), linTrack(std::move(lTrack)) {}

  BilloirTrack(const BilloirTrack& arg) = default;

  const input_track_t* originalTrack;
  Acts::LinearizedTrack linTrack;
  double chi2;
  Jacobian DiMat;                               // position jacobian
  Acts::ActsMatrix<Acts::eBoundSize, 3> EiMat;  // momentum jacobian
  Acts::ActsSymMatrix<3> CiMat;                 //  = EtWmat * Emat (see below)
  Acts::ActsMatrix<4, 3> BiMat;                 //  = DiMat^T * Wi * EiMat
  Acts::ActsSymMatrix<3> CiInv;                 //  = (EiMat^T * Wi * EiMat)^-1
  Acts::Vector3 UiVec;                          //  = EiMat^T * Wi * dqi
  Acts::ActsMatrix<4, 3> BCiMat;                //  = BiMat * Ci^-1
  Acts::BoundVector deltaQ;
};

/// @struct BilloirVertex
///
/// @brief Struct to cache vertex-specific matrix operations in Billoir fitter
struct BilloirVertex {
  // Amat  = sum{DiMat^T * Wi * DiMat}
  Acts::SymMatrix4 Amat = Acts::SymMatrix4::Zero();
  // Tvec  = sum{DiMat^T * Wi * dqi}
  Acts::Vector4 Tvec = Acts::Vector4::Zero();
  // BCBmat = sum{BiMat * Ci^-1 * BiMat^T}
  Acts::SymMatrix4 BCBmat = Acts::SymMatrix4::Zero();
  // BCUvec = sum{BiMat * Ci^-1 * UiVec}
  Acts::Vector4 BCUvec = Acts::Vector4::Zero();
};

}  // end anonymous namespace

template <typename input_track_t, typename linearizer_t>
Acts::Result<Acts::Vertex<input_track_t>>
Acts::FullBilloirVertexFitter<input_track_t, linearizer_t>::fit(
    const std::vector<const input_track_t*>& paramVector,
    const linearizer_t& linearizer,
    const VertexingOptions<input_track_t>& vertexingOptions,
    State& state) const {
  double chi2 = std::numeric_limits<double>::max();
  double newChi2 = 0;
  unsigned int nTracks = paramVector.size();

  if (nTracks == 0) {
    return Vertex<input_track_t>(Vector3(0., 0., 0.));
  }
  
  //if (nTracks != 2) {
  //  return Vertex<input_track_t>(Vector3(0., 0., 0.));
  //}

  // Set number of degrees of freedom
  // ndf = (5-3) * nTracks - 3;
  int ndf = 2 * nTracks - 3;
  if (nTracks < 2) {
    ndf = 1;
  }

  // Determine if we do contraint fit or not by checking if an
  // invertible non-zero constraint vertex covariance is given
  bool isConstraintFit = false;
  if (vertexingOptions.vertexConstraint.covariance().determinant() != 0) {
    isConstraintFit = true;
    ndf += 3;
  }

  std::vector<BilloirTrack<input_track_t>> billoirTracks;

  std::vector<Vector3> trackMomenta;

  Vector4 linPoint(vertexingOptions.vertexConstraint.fullPosition());

  Vertex<input_track_t> fittedVertex;

  for (int nIter = 0; nIter < m_cfg.maxIterations; ++nIter) {
    billoirTracks.clear();

    newChi2 = 0;

    BilloirVertex billoirVertex;
    int iTrack = 0;
    // iterate over all tracks
    
    std::cout << "loop over input tracks -----------------------" << std::endl;
    
    for (const input_track_t* trackContainer : paramVector) {
      const auto& trackParams = extractParameters(*trackContainer);
      if (nIter == 0) {
        double phi = trackParams.parameters()[BoundIndices::eBoundPhi];
        double theta = trackParams.parameters()[BoundIndices::eBoundTheta];
        double qop = trackParams.parameters()[BoundIndices::eBoundQOverP];
        
        //if (qop > 0) qop = 0.010000000000;
        //else qop = -0.010000000000;
        
        trackMomenta.push_back(Vector3(phi, theta, qop));
               
        //std::cout << "iter 0, qop " << qop << std::endl;
      }

      auto result = linearizer.linearizeTrack(
          trackParams, linPoint, vertexingOptions.geoContext,
          vertexingOptions.magFieldContext, state.linearizerState);
      if (result.ok()) {
        const auto& linTrack = *result;
        const auto& parametersAtPCA = linTrack.parametersAtPCA;
        double d0 = parametersAtPCA[BoundIndices::eBoundLoc0];
        double z0 = parametersAtPCA[BoundIndices::eBoundLoc1];
        double phi = parametersAtPCA[BoundIndices::eBoundPhi];
        double theta = parametersAtPCA[BoundIndices::eBoundTheta];
        double qOverP = parametersAtPCA[BoundIndices::eBoundQOverP];
        //if (qOverP > 0) qOverP = 0.010000000000;
        //else qOverP = -0.010000000000;
        
        // calculate f(V_0,p_0)  f_d0 = f_z0 = 0
        double fPhi = trackMomenta[iTrack][0];
        double fTheta = trackMomenta[iTrack][1];
        double fQOvP =  trackMomenta[iTrack][2];
        BilloirTrack<input_track_t> currentBilloirTrack(trackContainer,
                                                        linTrack);
                                                        
        //std::cout << "if result of linearizer is ok, qop " << qOverP << " fqop "<< fQOvP << std::endl;                                                

        currentBilloirTrack.deltaQ << d0, z0, phi - fPhi, theta - fTheta,
            qOverP - fQOvP, 0;

        // position jacobian (D matrix)
        ActsMatrix<eBoundSize, 4> Dmat;
        Dmat = linTrack.positionJacobian;

        // momentum jacobian (E matrix)
        ActsMatrix<eBoundSize, 3> Emat;
        Emat = linTrack.momentumJacobian;
        // cache some matrix multiplications
        ActsMatrix<4, eBoundSize> DtWmat;
        ActsMatrix<3, eBoundSize> EtWmat;
        BoundSymMatrix Wi = linTrack.weightAtPCA;

        DtWmat = Dmat.transpose() * Wi;
        EtWmat = Emat.transpose() * Wi;

        // compute billoir tracks
        currentBilloirTrack.DiMat = Dmat;
        currentBilloirTrack.EiMat = Emat;
        currentBilloirTrack.CiMat = EtWmat * Emat;
        currentBilloirTrack.BiMat = DtWmat * Emat;  // DiMat^T * Wi * EiMat
        currentBilloirTrack.UiVec =
            EtWmat * currentBilloirTrack.deltaQ;  // EiMat^T * Wi * dqi
        currentBilloirTrack.CiInv =
            (EtWmat * Emat).inverse();  // (EiMat^T * Wi * EiMat)^-1


        auto test_mat = (EtWmat * Emat); //;.inverse();
        
        std::cout << "Wi +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
        std::cout << Wi  << std::endl;
        
        std::cout << "Emat +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
        std::cout << Emat  << std::endl;
        
        std::cout << "CiInv determinante oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo " << test_mat.determinant()  << std::endl;
        std::cout << "CiInv inverse +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
        std::cout << test_mat.inverse()  << std::endl;

        // sum up over all tracks
        billoirVertex.Tvec +=
            DtWmat * currentBilloirTrack.deltaQ;  // sum{DiMat^T * Wi * dqi}
        billoirVertex.Amat += DtWmat * Dmat;      // sum{DiMat^T * Wi * DiMat}

        // remember those results for all tracks
        currentBilloirTrack.BCiMat =
            currentBilloirTrack.BiMat *
            currentBilloirTrack.CiInv;  // BCi = BiMat * Ci^-1

        // and some summed results
        billoirVertex.BCUvec +=
            currentBilloirTrack.BCiMat *
            currentBilloirTrack.UiVec;  // sum{BiMat * Ci^-1 * UiVec}
        billoirVertex.BCBmat +=
            currentBilloirTrack.BCiMat *
            currentBilloirTrack.BiMat
                .transpose();  // sum{BiMat * Ci^-1 * BiMat^T}

        billoirTracks.push_back(currentBilloirTrack);
        ++iTrack;
      } else {
        return result.error();
      }
    }  // end loop tracks

    // calculate delta (billoirFrameOrigin-position), might be changed by the
    // beam-const
    // Vdel = Tvec-sum{BiMat*Ci^-1*UiVec}
    Vector4 Vdel = billoirVertex.Tvec - billoirVertex.BCUvec;
    
    SymMatrix4 VwgtMat =
        billoirVertex.Amat -
        billoirVertex.BCBmat;  // VwgtMat = Amat-sum{BiMat*Ci^-1*BiMat^T}
        
    if (isConstraintFit) {
      // this will be 0 for first iteration but != 0 from second on
      Vector4 posInBilloirFrame =
          vertexingOptions.vertexConstraint.fullPosition() - linPoint;

      Vdel += vertexingOptions.vertexConstraint.fullCovariance().inverse() *
              posInBilloirFrame;
      VwgtMat += vertexingOptions.vertexConstraint.fullCovariance().inverse();
    }

    // cov(deltaV) = VwgtMat^-1
    
    std::cout << "VwgtMat ------ " << std::endl;
    std::cout << VwgtMat << std::endl;
    
    /*std::vector<std::vector<double>> matrix(4, std::vector<double> (4));
    /matrix = {
        {VwgtMat(0,0),VwgtMat(0,1),VwgtMat(0,2), VwgtMat(0,3)},
        {VwgtMat(1,0),VwgtMat(1,1),VwgtMat(1,2), VwgtMat(1,3)},
        {VwgtMat(2,0),VwgtMat(2,1),VwgtMat(2,2), VwgtMat(2,3)},
        {VwgtMat(3,0),VwgtMat(3,1),VwgtMat(3,2), VwgtMat(3,3)}
    };*/

    //printMatrixBernie(getInverseBernie(matrix));
    
    std::cout << "determinante xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << VwgtMat.determinant()  << std::endl;
    //std::cout << VwgtMat.isInvertible() << std::endl;
    Eigen::FullPivLU<SymMatrix4> lu(VwgtMat);
    std::cout << lu.isInvertible() << std::endl;
    
    SymMatrix4 covDeltaVmat = VwgtMat.inverse();
    
    std::cout << "covDeltaVmat ------ " << std::endl;
    std::cout << covDeltaVmat << std::endl;
    
    
    
    // deltaV = cov_(deltaV) * Vdel;
    Vector4 deltaV = covDeltaVmat * Vdel;
    //--------------------------------------------------------------------------------------
    // start momentum related calculations

    std::vector<std::optional<BoundSymMatrix>> covDeltaPmat(nTracks);

    iTrack = 0;
    
    std::cout << "loop over billoirt tracks -----------------------" << std::endl;
    
    for (auto& bTrack : billoirTracks) {
      Vector3 deltaP =
          (bTrack.CiInv) * (bTrack.UiVec - bTrack.BiMat.transpose() * deltaV);

      // update track momenta
      trackMomenta[iTrack] += deltaP;
      
      //std::cout << "track momenta " << trackMomenta[iTrack] << std::endl;

      // correct for 2PI / PI periodicity
      const auto correctedPhiTheta = detail::normalizePhiTheta(
          trackMomenta[iTrack][0], trackMomenta[iTrack][1]);
      trackMomenta[iTrack][0] = correctedPhiTheta.first;
      trackMomenta[iTrack][1] = correctedPhiTheta.second;
      
      //std::cout << "track momenta corr " << trackMomenta[iTrack] << std::endl;

      // calculate 5x5 covdelta_P matrix
      // d(d0,z0,phi,theta,qOverP, t)/d(x,y,z,phi,theta,qOverP,
      // t)-transformation matrix
      ActsMatrix<eBoundSize, 7> transMat;
      transMat.setZero();
      transMat(0, 0) = bTrack.DiMat(0, 0);
      transMat(0, 1) = bTrack.DiMat(0, 1);
      transMat(1, 0) = bTrack.DiMat(1, 0);
      transMat(1, 1) = bTrack.DiMat(1, 1);
      transMat(1, 2) = 1.;
      transMat(2, 3) = 1.;
      transMat(3, 4) = 1.;
      transMat(4, 5) = 1.;
      transMat(5, 6) = 1.;

      // some intermediate calculations to get 5x5 matrix
      // cov(V,V), 4x4 matrix
      SymMatrix4 VVmat = covDeltaVmat;

      // cov(V,P)
      ActsMatrix<4, 3> VPmat = bTrack.BiMat;

      // cov(P,P), 3x3 matrix
      ActsSymMatrix<3> PPmat;
      PPmat = bTrack.CiInv +
              bTrack.BCiMat.transpose() * covDeltaVmat * bTrack.BCiMat;

      ActsSymMatrix<7> covMat;
      covMat.setZero();
      covMat.block<4, 4>(0, 0) = VVmat;
      covMat.block<4, 3>(0, 4) = VPmat;
      covMat.block<3, 4>(4, 0) = VPmat.transpose();

      covMat.block<3, 3>(4, 4) = PPmat;

      // covdelta_P calculation
      covDeltaPmat[iTrack] = transMat * covMat * transMat.transpose();
      // Calculate chi2 per track.
      bTrack.chi2 =
          ((bTrack.deltaQ - bTrack.DiMat * deltaV - bTrack.EiMat * deltaP)
               .transpose())
              .dot(bTrack.linTrack.weightAtPCA *
                   (bTrack.deltaQ - bTrack.DiMat * deltaV -
                    bTrack.EiMat * deltaP));
      newChi2 += bTrack.chi2;

      ++iTrack;
    }

    if (isConstraintFit) {
      // last term will also be 0 again but only in the first iteration
      // = calc. vtx in billoir frame - (    isConstraintFit pos. in billoir
      // frame )

      Vector4 deltaTrk =
          deltaV -
          (vertexingOptions.vertexConstraint.fullPosition() - linPoint);

      newChi2 +=
          (deltaTrk.transpose())
              .dot(
                  vertexingOptions.vertexConstraint.fullCovariance().inverse() *
                  deltaTrk);
    }

    if (!std::isnormal(newChi2)) {
      return VertexingError::NumericFailure;
    }

    // assign new linearization point (= new vertex position in global frame)
    linPoint += deltaV;
    
    //std::cout << "old chi2 " << chi2 << ", new chi2 " << newChi2 << std::endl;
    
    std::cout << "lin point ---------------------- " << std::endl;
    std::cout << linPoint << std::endl;
    
    //std::cout << "covariance matrix "<< covDeltaVmat.diagonal()[0] << std::endl;
    
    if (newChi2 < chi2) {
      chi2 = newChi2;
      
      //std::cout << " adding this track!" << std::endl;

      Vector4 vertexPos(linPoint);

      fittedVertex.setFullPosition(vertexPos);
      
      
      
      
      fittedVertex.setFullCovariance(covDeltaVmat);
      fittedVertex.setFitQuality(chi2, ndf);

      std::vector<TrackAtVertex<input_track_t>> tracksAtVertex;

      std::shared_ptr<PerigeeSurface> perigee =
          Surface::makeShared<PerigeeSurface>(
              VectorHelpers::position(vertexPos));

      iTrack = 0;
      for (auto& bTrack : billoirTracks) {
        // new refitted trackparameters
        BoundVector paramVec = BoundVector::Zero();
        paramVec[eBoundPhi] = trackMomenta[iTrack](0);
        paramVec[eBoundTheta] = trackMomenta[iTrack](1);
        paramVec[eBoundQOverP] = trackMomenta[iTrack](2);
        BoundTrackParameters refittedParams(perigee, paramVec,
                                            covDeltaPmat[iTrack]);
        TrackAtVertex<input_track_t> trackVx(bTrack.chi2, refittedParams,
                                             bTrack.originalTrack);
        tracksAtVertex.push_back(std::move(trackVx));
        ++iTrack;
      }
      fittedVertex.setTracksAtVertex(tracksAtVertex);
    }
    else {
    //std::cout << "chi2 is larger than before!" << std::endl;
    }
  }  // end loop iterations
  return fittedVertex;
}
