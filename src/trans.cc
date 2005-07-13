// ----------------------------------------------------------------------------
// Title: trans
// ----------------------------------------------------------------------------
// Author: Matthias Wangler
//         mw@imbi.uni-freiburg.de
//         Institute of Med. Biometry and Med. Computer Science
//         Stefan-Meier-Strasse 26, D-79104 Freiburg,
//         http://www.imbi.uni-freiburg.de
// ----------------------------------------------------------------------------
// Description: Computes the transition matrices 
// ----------------------------------------------------------------------------
// Required: STL (C++ Standard-Template-Library)
// ----------------------------------------------------------------------------
// Usage in an R-Program:
//
// .C("trans", as.integer(states), 
//             as.integer(length(state.names)),
//             as.integer(model$jumps[,1]), 
//             as.integer(model$jumps[,2]), 
//             as.integer(nrow(model$jumps)),
//             nj = as.integer(nrjumps[,3]), 
//             nr = as.integer(nr.before),
//             as.double(times), 
//             as.integer(length(times)),
//             as.integer(match(as.character(observ$from),state.names)),
//             as.integer(match(as.character(observ$to),state.names)),
//             as.double(observ$time), 
//             as.integer(length(observ$from)),
//             ma = as.double(matrices), 
//             as.integer(nrow(model$trans)), 
//             as.integer(ncol(model$trans)), 
//             as.integer(length(times)) )
// ----------------------------------------------------------------------------
// Value:
// 
// nj: an vector of integers which is holding the total number of jumnps for
//     each possible transition
// ma: an double vector which is holding the transition matrices        
// ----------------------------------------------------------------------------
// Notes: -
// ----------------------------------------------------------------------------
// Example: see trans.R
// ----------------------------------------------------------------------------
// License: GPL 2
//-----------------------------------------------------------------------------
// History: 30.08.2004, Matthias Wangler
//                      the first version
//          16.02.2005, Matthias Wangler
//                      replace long with int
// ----------------------------------------------------------------------------

#include <math.h>
extern "C" {
  void trans( int* states,
	      int* nr_states, 
	      int* from, 
	      int* to, 
	      int* nr_from_to, 
	      int* nr_jumps,
	      int* nr_before,           
	      double* times, 
	      int* nr_times,
	      int* observ_from, 
	      int* observ_to, 
	      double* observ_time,
	      int* observ_time_point, 
	      int* nr_observ,
	      double* a, 
	      int* nrow, 
	      int* ncol ) {

    int  CensState = states[*nr_states-1];  // the censoring state

    int no   = *nr_observ;
    int nft  = *nr_from_to;
    int rows = *nrow;
    int cols = *ncol;
    int nt   = *nr_times;

    // compute the number of jumps and store it in the array of the transition matrices
    for(int i=0 ; i < no; ++i) {      // loop over the observations (jumps)
      // observ_from[i]: state from where the observed jump occurs   
      // observ_to[i]: state to which the observed jump occurs
      // observ_time_point[i]: the observed time point
	
      int k = observ_time_point[i]-1;        
      int r = observ_from[i]-1;        
      int c = observ_to[i]-1;

      // add the jump to the total number of jumps from state 'observ_from[i]' to state 'observ_to[i]'
      for( int j = 0; j < nft; ++j ) {
	if( observ_from[i] == from[j] && observ_to[i] == to[j] ) {
	  nr_jumps[j] += 1;
	}
      }

      if( observ_to[i] != CensState ) {   // jumps to the censoring state are ignored
	// add this jump in the k'th matrix at column c and row r	             
	a[(r+(c*rows)) + k*(rows*cols)] += 1;
      }

      if( observ_time_point[i] < nt ) {      
	nr_before[(k+1) + (r*nt)] -= 1;
	nr_before[(k+1) + (c*nt)] += 1; 
      }
    }

    for(int t=0; t < nt; ++t) {     
      for( int j = 0; j < rows; ++j ) { 
	if( t > 0 ) {
	  nr_before[t + (j*nt)] = nr_before[(t-1) + (j*nt)] +  nr_before[t + (j*nt)];
	}

	if(nr_before[t + (nt*j)] > 0 ) {	   
	  for( int c = 0; c < cols; ++c ) {             	     
	    a[(j+(c*rows)) + t*(rows*cols)] = a[(j+(c*rows)) + t*(rows*cols)] / nr_before[t + (j*nt)];	   
	  }
	}
      }

      // compute the diagonal elements: the sum over each row must be 1
      for( int r = 0; r < rows; ++r ) {
        double sumrow = 0;

        for( int c = 0; c < cols; ++c ) {
	  if( c != r ) {
	    sumrow += a[(r+(c*rows)) + t*(rows*cols)];
	  }	 
        }

	a[(r+(r*rows)) + t*(rows*cols)] = (double)(1 - sumrow); 
      }    
    }
  }
}

