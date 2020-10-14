// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
    
IntegerVector which_eq(IntegerVector x, int j) {
  // Returns a vector of the same length with non-zero elements indicating which elements of x = j.
  int n = x.size();
  IntegerVector check(n);
  for (int i = 0; i < n; i++) {
    if (x[i]!=j) {
      check[i] = 0;
    } else {
      check[i] = i+1;
    }
  }
  return check;
}

IntegerVector which_ge(IntegerVector x, int j) {
  // Returns a vector of the same length with non-zero elements indicating which elements of x >= j.
  int n = x.size();
  IntegerVector check(n);
  for (int i = 0; i < n; i++) {
    if (x[i]<j) {
      check[i] = 0;
    } else {
      check[i] = i+1;
    }
  }
  return check;
}

rowvec which_ge_arma(rowvec x, int j) {
  // Returns a vector of the same length with non-zero elements indicating which elements of x >= j.
  int n = x.size();
  rowvec check(n);
  for (int i = 0; i < n; i++) {
    if (x[i]<j) {
      check[i] = 0;
    } else {
      check[i] = i+1;
    }
  }
  return check;
}

double prod_Cpp(NumericVector x) {
  // Returns the product of all elements in the vector x.
  int n = x.size();
  double x_prod = x[0];
  for (int i = 1; i < n; i++) {
    x_prod *= x[i];
  }
  return x_prod;
}

rowvec mvrnorm_Cpp(int n, vec mu, mat Sigma) {
  int ncols = Sigma.n_cols;
  rowvec Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(Sigma);
}

mat mnl_prob_Cpp(colvec beta, int y, rowvec incl_alts, mat X) {
// mat mnl_prob_Cpp(colvec beta, int y, IntegerVector incl_alts, mat X) {
  // Returns the closed-form MNL choice probability for the given choice task.
  //   beta = coefficients.
  //   y = Picked alternative.
  //   incl_alts = Vector of included alternatives.
  //   X = Current X.
  
  // Count the number of included alternatives.
  int l_incl_alts = incl_alts.size();
  int n_incl_alts = 0;
  for (int i = 0; i < l_incl_alts; i++) {
    if (incl_alts[i]!=0) {
      n_incl_alts += 1;
    }
  }
  
  // Create a vector of row indices for the included alternatives.
  uvec incl_alts_sub(n_incl_alts);
  int j = 0;
  for (int i = 0; i < l_incl_alts; i++) {
    if (incl_alts[i]!=0) {
      incl_alts_sub[j] = incl_alts[i]-1;
      j += 1;
    }
  }
  
  // Create X_s using X for the included alternatives.
  mat X_s = X.rows(incl_alts_sub);
  
  // Compute the closed-form MNL choice probability.
  rowvec x_pick = X.row(y-1);
  mat Numerator = exp(x_pick*beta);
  mat Denominator = zeros<mat>(1,1);
  for (int j = 0; j < n_incl_alts; j++) {
    rowvec x_j = X_s.row(j);
    Denominator += exp(x_j*beta);
  }
  
  // Return the choice probability.
  return (Numerator/Denominator);
}

double Xbeta_mnl_prob_Cpp(int y, IntegerVector incl_alts, NumericVector Xbeta) {
  // Returns the closed-form MNL choice probability for the given choice task.
  //   y = Picked alternative.
  //   incl_alts = Vector of included alternatives.
  //   Xbeta = Current X%*%beta.
  
  // Count the number of included alternatives.
  int l_incl_alts = incl_alts.size();
  int n_incl_alts = 0;
  for (int i = 0; i < l_incl_alts; i++) {
    if (incl_alts[i]==0) {
      n_incl_alts += 0;
    } else {
      n_incl_alts += 1;
    }
  }
  
  // Create Xbeta_s using Xbeta for the included alternatives.
  NumericVector Xbeta_s(n_incl_alts);
  IntegerVector incl_alts_s(n_incl_alts);
  int j = 0;
  for (int i = 0; i < l_incl_alts; i++) {
    if (incl_alts[i]==0) {
      j += 0;
    } else {
      Xbeta_s[j] = Xbeta[i];
      incl_alts_s[j] = incl_alts[i];
      j += 1;
    }
  }
  
  // Identify the pick in the new subset.
  int y_s = sum(which_eq(incl_alts_s,y));
  
  // Compute the closed-form MNL choice probability.
  double Numerator = exp(Xbeta_s[y_s-1]);
  double Denominator = 0;
  for (int i = 0; i < n_incl_alts; i++) {
    Denominator += exp(Xbeta_s[i]);
  }
  
  // Return the choice probability.
  return (Numerator/Denominator);
}

// [[Rcpp::export]]
List conj_llmnl_Cpp(int nscns, int nalts, int nlvls, int nvars, IntegerVector y, mat X_ll, IntegerMatrix X_rel, 
                    IntegerMatrix S, NumericVector oldtaus, colvec betad, colvec betac, int out_ind) {
  // Evaluate the likelihood for old and canidate beta draws.
  //   nscns = Number of choice tasks.
  //   nalts = Number of alternatives in each choice task.
  //   nlvls = Number of levels for non-brand attributes.
  //   nvars = Number of attribute levels (excluding reference levels).
  //   y = Vector of choices.
  //   X_ll = Design matrix.
  //   X_rel = Design matrix.
  //   S = Levels matrix.
  //   oldtaus = Current draws of oldtaus.
  //   betad = Old values of beta.
  //   betac = Candidates value of beta.
  //   out_ind = Indicates an outside good.
  
  // Initialize a list and matrices for the output: logold and lognew.
  List out(3);
  mat logold = zeros<mat>(1,1);
  mat lognew = zeros<mat>(1,1);
  
  // Likelihood scenario loop.
  for (int scn = 0; scn < nscns; scn++) {
    // Indicate attribute relevance for each alternative (outside option, when present, is always attribute relevant).
    rowvec attr_rel_scn(nalts);
    if (out_ind==0) {
      IntegerMatrix S_scn = S(Range(((scn+1)*nalts-nalts),((scn+1)*nalts-1)),_);
      for (int alt = 0; alt < nalts; alt++) {
        NumericVector attr_rel_alt(nlvls);
        for (int lvl = 0; lvl < nlvls; lvl++) {
          attr_rel_alt[lvl] = 1 - oldtaus[lvl]*(S_scn(alt,lvl));
        }
        attr_rel_scn[alt] = prod_Cpp(attr_rel_alt);
      }
    } else {
      IntegerMatrix S_scn = S(Range(((scn+1)*(nalts-1)-(nalts-1)),((scn+1)*(nalts-1)-1)),_);
      for (int alt = 0; alt < (nalts-1); alt++) {
        NumericVector attr_rel_alt(nlvls);
        for (int lvl = 0; lvl < nlvls; lvl++) {
          attr_rel_alt[lvl] = 1 - oldtaus[lvl]*(S_scn(alt,lvl));
        }
        attr_rel_scn[alt] = prod_Cpp(attr_rel_alt);
      }
      attr_rel_scn[nalts-1] = 1;
    }
    
    // Log likelihood with old and candidate beta draws.
    rowvec incl_alts = which_ge_arma(attr_rel_scn,1);
    logold += log(mnl_prob_Cpp(betad,y[scn],incl_alts,X_ll.submat((scn+1)*nalts-nalts,0,(scn+1)*nalts-1,nvars-1)));
    lognew += log(mnl_prob_Cpp(betac,y[scn],incl_alts,X_ll.submat((scn+1)*nalts-nalts,0,(scn+1)*nalts-1,nvars-1)));
  }
  
  // Return logold and lognew.
  out[0] = logold;
  out[1] = lognew;
  return out;
}

// [[Rcpp::export]]
NumericVector draw_tau_Cpp(int resp, int nresp, int nscns, int nalts, int nlvls, IntegerVector y, IntegerMatrix S, 
                           NumericVector Xbeta, NumericVector oldtaus, NumericMatrix oldtheta, int out_ind, int het_ind) {
  // Draw tau, the screening parameters.
  //   resp = Current respondent.
  //   nresp = Number of respondents.
  //   nscns = Number of choice tasks.
  //   nalts = Number of alternatives in each choice task.
  //   nlvls = Number of levels for attributes.
  //   y = Vector of choices.
  //   S = Levels matrix.
  //   Xbeta = X%*%beta.
  //   oldtaus = Current draws of oldtaus to update.
  //   oldtheta = Current draws of oldthetas.
  //   out_ind = Indicates an outside good.
  //   het_ind = Indicates heterogeneity in Theta.
  
  // Clone the draws to update.
  NumericVector oldtaus_new = clone(oldtaus);
  
  // Tau non-brand attribute-levels loop.
  for (int lvl = 0; lvl < nlvls; lvl++) {
    // Counting picks including lvl scenario loop.
    IntegerMatrix incl_ind(nscns,1);
    for (int scn = 0; scn < nscns; scn++) {
      if (out_ind==0) {
        // Levels matrix.
        IntegerMatrix S_scn = S(Range(((scn+1)*nalts-nalts),((scn+1)*nalts-1)),_);
        
        // Was lvl in the picked alt?
        if ( sum(which_eq(which_eq(S_scn(_,lvl),1),y[scn]))==0 ) {
          incl_ind(scn,0) = 0;
        } else {
          incl_ind(scn,0) = 1;
        }
      } else {
        // Levels matrix.
        IntegerMatrix S_scn = S(Range(((scn+1)*(nalts-1)-(nalts-1)),((scn+1)*(nalts-1)-1)),_);
        
        // Was lvl in the picked alt?
        if ( sum(which_eq(which_eq(S_scn(_,lvl),1),y[scn]))==0 ) {
          incl_ind(scn,0) = 0;
        } else {
          incl_ind(scn,0) = 1;
        }
      }
    }
    
    // Set or draw tau based on if the picks were brand irrelevant and if the level was in the picks.
    if ( sum((incl_ind(_,0)==1))>=1 ) {
      oldtaus_new[lvl] = 0;
    } else {
      // Tau zero and one draws.
      NumericVector tau_zero = clone(oldtaus_new);
      NumericVector tau_one = clone(oldtaus_new);
      tau_zero[lvl] = 0;
      tau_one[lvl] = 1;
      
      // Drawing tau when lvl is included in picks scenario loop.
      double log_tau_zero = 0;
      double log_tau_one = 0;
      for (int scn = 0; scn < nscns; scn++) {
        if (out_ind==0) {
          // Level matrix.
          IntegerMatrix S_scn = S(Range(((scn+1)*nalts-nalts),((scn+1)*nalts-1)),_);
          
          // Indicate attribute relevance for each alternative.
          IntegerVector attr_rel_scn_zero(nalts);
          IntegerVector attr_rel_scn_one(nalts);
          for (int alt = 0; alt < nalts; alt++) {
            NumericVector attr_rel_alt_zero(nlvls);
            NumericVector attr_rel_alt_one(nlvls);
            for (int lvl_2 = 0; lvl_2 < nlvls; lvl_2++) {
              attr_rel_alt_zero[lvl_2] = 1 - tau_zero[lvl_2]*(S_scn(alt,lvl_2));
              attr_rel_alt_one[lvl_2] = 1 - tau_one[lvl_2]*(S_scn(alt,lvl_2));
            }
            attr_rel_scn_zero[alt] = prod_Cpp(attr_rel_alt_zero);
            attr_rel_scn_one[alt] = prod_Cpp(attr_rel_alt_one);
          }
          
          // Xbeta for current scenario.
          NumericVector Xbeta_scn = Xbeta[Range(((scn+1)*nalts-nalts),((scn+1)*nalts-1))];
          
          // Log likelihood with old and candidate tau draws.
          IntegerVector incl_alts_zero = which_ge(attr_rel_scn_zero,1);
          IntegerVector incl_alts_one = which_ge(attr_rel_scn_one,1);
          if (sum(which_eq(incl_alts_zero,y[scn]))==0) {
            log_tau_zero += 0;
          } else {
            log_tau_zero += log(Xbeta_mnl_prob_Cpp(sum(which_eq(incl_alts_zero,y[scn])),incl_alts_zero,Xbeta_scn));
          }
          if (sum(which_eq(incl_alts_one,y[scn]))==0) {
            log_tau_one += 0;
          } else {
            log_tau_one += log(Xbeta_mnl_prob_Cpp(sum(which_eq(incl_alts_one,y[scn])),incl_alts_one,Xbeta_scn));
          }
        } else {
          // Level matrix.
          IntegerMatrix S_scn = S(Range(((scn+1)*(nalts-1)-(nalts-1)),((scn+1)*(nalts-1)-1)),_);
          
          // Indicate attribute relevance for each alternative (outside option always attribute relevant).
          IntegerVector attr_rel_scn_zero(nalts);
          IntegerVector attr_rel_scn_one(nalts);
          for (int alt = 0; alt < (nalts-1); alt++) {
            NumericVector attr_rel_alt_zero(nlvls);
            NumericVector attr_rel_alt_one(nlvls);
            for (int lvl_2 = 0; lvl_2 < nlvls; lvl_2++) {
              attr_rel_alt_zero[lvl_2] = 1 - tau_zero[lvl_2]*(S_scn(alt,lvl_2));
              attr_rel_alt_one[lvl_2] = 1 - tau_one[lvl_2]*(S_scn(alt,lvl_2));
            }
            attr_rel_scn_zero[alt] = prod_Cpp(attr_rel_alt_zero);
            attr_rel_scn_one[alt] = prod_Cpp(attr_rel_alt_one);
          }
          attr_rel_scn_zero[nalts-1] = 1;
          attr_rel_scn_one[nalts-1] = 1;
          
          // Xbeta for current scenario.
          NumericVector Xbeta_scn = Xbeta[Range(((scn+1)*nalts-nalts),((scn+1)*nalts-1))];
          
          // Log likelihood with old and candidate tau draws.
          IntegerVector incl_alts_zero = which_ge(attr_rel_scn_zero,1);
          IntegerVector incl_alts_one = which_ge(attr_rel_scn_one,1);
          if (sum(which_eq(incl_alts_zero,y[scn]))==0) {
            log_tau_zero += 0;
          } else {
            log_tau_zero += log(Xbeta_mnl_prob_Cpp(sum(which_eq(incl_alts_zero,y[scn])),incl_alts_zero,Xbeta_scn));
          }
          if (sum(which_eq(incl_alts_one,y[scn]))==0) {
            log_tau_one += 0;
          } else {
            log_tau_one += log(Xbeta_mnl_prob_Cpp(sum(which_eq(incl_alts_one,y[scn])),incl_alts_one,Xbeta_scn));
          }
        }
      }
      
      // Compute screening probability.
      double var_theta = 0;
      if (het_ind==0) {
        var_theta = ((exp(log_tau_one)*oldtheta[lvl])/(exp(log_tau_one)*oldtheta[lvl]+exp(log_tau_zero)*(1-oldtheta[lvl])));
      } else {
        var_theta = ((exp(log_tau_one)*oldtheta(resp-1,lvl))/(exp(log_tau_one)*oldtheta(resp-1,lvl)+exp(log_tau_zero)*(1-oldtheta(resp-1,lvl))));
      }
      
      // Draw tau.
      oldtaus_new[lvl] = R::rbinom(1,var_theta);
    }
  }
  return oldtaus_new;
}

// [[Rcpp::export]]
List draw_omega_theta_Cpp(int nresp, int nlvls, int nscov, mat Z, mat D, mat oldOmega, double ostep, int onaccept, 
                          mat omegabar, mat Vomegabar, mat Vomegabari, IntegerMatrix oldtaus, mat oldtheta) {
  // Draw Omega and theta using homogeneous logits.
  //   nresp = Number of respondents.
  //   nlvls = Number of levels for non-brand attributes.
  //   nscov = Number of screening covariates.
  //   Z = Covariates for the upper-level screening model.
  //   D = "Design matrix" for the upper-level screening model.
  //   oldOmega = Current draws of oldOmega to update to update.
  //   ostep = Current step size for the Omega draw.
  //   onaccept = Number of times Omega has been accepted.
  //   omegabar = Means for normal prior on Omega.
  //   Vomegabar = Covariance matrix for normal prior on Omega.
  //   Vomegabari = chol2inv(chol(Vomegabar)).
  //   oldtaus = Current draws of oldtaus.
  //   oldtheta = Current draws of oldtheta to update.
  
  // Clone the draws to update.
  mat oldOmega_new = conv_to<mat>::from(oldOmega);
  int onaccept_new = onaccept;
  mat oldtheta_new = conv_to<mat>::from(oldtheta);
  
  // Initialize a list for the output: oldOmega_new, onaccept_new, and oldtheta_new.
  List out(3);
  
  // Omega and theta non-brand attribute-levels loop.
  for (int lvl = 0; lvl < nlvls; lvl++) {
    // Omega old and candidate draws.
    rowvec omegad_r = oldOmega_new.row(lvl);
    rowvec omegac_r = omegad_r + mvrnorm_Cpp(1,zeros<vec>(nscov+1),(ostep*Vomegabar));
    
    // Log likelihood for old and candidate draws.
    mat logoldlike = zeros<mat>(1,1);
    mat lognewlike = zeros<mat>(1,1);
    colvec omegad_c = omegad_r.t();
    colvec omegac_c = omegac_r.t();
    rowvec incl_alts = zeros<rowvec>(2);
    incl_alts[0] = 1;
    incl_alts[1] = 2;
    for (int resp = 0; resp < nresp; resp++) {
      logoldlike += log(mnl_prob_Cpp(omegad_c,oldtaus(resp,lvl)+1,incl_alts,D.submat((resp+1)*2-2,0,(resp+1)*2-1,nscov)));
      lognewlike += log(mnl_prob_Cpp(omegac_c,oldtaus(resp,lvl)+1,incl_alts,D.submat((resp+1)*2-2,0,(resp+1)*2-1,nscov)));
    }
    
    // Log priors on Omega.
    mat omegabar_Cpp = conv_to<mat>::from(omegabar);
    rowvec omegabar_r = omegabar_Cpp.row(lvl);
    colvec omegabar_c = omegabar_Cpp.row(lvl).t();
    mat logoldprior = -0.5*(omegad_r-omegabar_r)*Vomegabari*(omegad_c-omegabar_c);
    mat lognewprior = -0.5*(omegac_r-omegabar_r)*Vomegabari*(omegac_c-omegabar_c);
    
    // Omega log posteriors.
    mat logoldpost = logoldlike + logoldprior;
    mat lognewpost = lognewlike + lognewprior;
    
    // Compare the old and candidate posteriors and compute kappa (there is no second-stage prior).
    mat diff = exp((lognewpost) - (logoldpost));
    mat kappa(1,1);
    if (diff.has_nan() || diff.has_inf()) {
      kappa = -1; // If the number doesn't exist, always reject.
    } else {
      kappa = min(ones<mat>(1,1),diff);
    }
    mat unif = randu<mat>(1,1);
    if (as_scalar(unif) < as_scalar(kappa)) {
      oldOmega_new.row(lvl) = omegac_r;
      onaccept_new += 1;
    }
    
    // Draw theta for lvl.
    rowvec oldOmega_r = oldOmega_new.row(lvl);
    for (int resp = 0; resp < nresp; resp++) {
      colvec z = Z.row(resp).t();
      oldtheta_new(resp,lvl) = 1 - as_scalar(exp(oldOmega_r*z)/(exp(oldOmega_r*z)+1));
    }
  }
  
  // Return oldOmega_new, onaccept_new, and oldtheta_new.
  out[0] = oldOmega_new;
  out[1] = onaccept_new;
  out[2] = oldtheta_new;
  return out;
}
