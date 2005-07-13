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
// Usage in an R-Program:
//
//  out <- .C("trans",
//             as.integer(length(states)-1),
//             as.integer(length(times)),
//             as.integer(length(observ.time)),
//             as.integer(model$transitions[,1]-1),
//             as.integer(model$transitions[,2]-1),
//             as.integer(nrow(model$transitions)),
//             nj = as.integer(nrtransitions[,3]),
//             as.double(observ.time),
//             as.integer(observ.event),
//             as.integer(observ.from),
//             as.integer(observ.to),
//             nr = as.integer(nr.before),
//             ma = as.double(matrices),
//             PACKAGE="changeLOS")
// ----------------------------------------------------------------------------
// Value:
// 
// nj: an vector of integers which is holding the total number of jumnps for
//     each possible transition
// nr: an vector to store the number in each state just before the transition 
//     times (risk)
// ma: an double vector which is holding the transition matrices        
// ----------------------------------------------------------------------------
// Notes: -
// ----------------------------------------------------------------------------
// Example: see trans.R
// ----------------------------------------------------------------------------
// License: GPL 2
//-----------------------------------------------------------------------------
#include <math.h>

using namespace std;

extern "C" {

  void trans(int* nr_states,
	     int* nr_time_points,
	     int* nr_observ,
	     int* from, 
	     int* to, 
	     int* nr_from_to, 
	     int* nr_jumps,
	     double* observ_time,
	     int* observ_event,
	     int* observ_from,
	     int* observ_to,               
	     double* risk,
	     double* trans) {
    
    int k = 0;
    int r = 0;
    int c = 0;
    int i = 0;
    int j = 0;
    int s = 0;
    
    int nft  = *nr_from_to;
    int no =  *nr_observ;
    int rows = *nr_states;
    int cols = *nr_states;
    
    int nt = *nr_time_points;
    
    double prev = observ_time[0];
  
  
    // compute the number of jumps and store it in the array of the transition matrices
    for(i=0; i < no; ++i) { // loop over the observations (jumps)
      r = observ_from[i];
      c = observ_to[i];
      
      if( prev != observ_time[i] && observ_event[i] == 1 ) {
	// increase the time point
	++k;
      }
      
      // add the jump to the total number of jumps from state 'observ_from[i]' to state 'observ_to[i]'
      for( s = 0; s < nft; ++s ) {
	if( r == from[s] && c == to[s] ) {
	  nr_jumps[s] += 1;
	}
      }
    
      if( observ_event[i] == 1 ) {	
	// add this jump in the k'th matrix at column c and row r	             	
	trans[(r+(c*rows)) + k*(rows*cols)] += 1;
      
      }
    
      if( k < (nt-1) ) {      
	risk[(k+1) + (r*nt)] -= 1;
	risk[(k+1) + (c*nt)] += 1; 
      }
    
      prev = observ_time[i];
    }
    
    for(k=0; k < nt; ++k) {          
      for( j = 0; j < rows; ++j ) { 	
	if( k > 0 ) {	  
	  risk[k + (j*nt)] = risk[(k-1) + (j*nt)] +  risk[k + (j*nt)];	
	}
	
	if(risk[k + (nt*j)] > 0 ) {	   	  
	  for( c = 0; c < cols; ++c ) {             	     	    
	    trans[(j+(c*rows)) + k*(rows*cols)] = trans[(j+(c*rows)) + k*(rows*cols)] / risk[k + (j*nt)];
	  }	
	}     
      }
      
      // compute the diagonal elements: the sum over each row must be 1 
      for( r = 0; r < rows; ++r ) {        
	double sumrow = 0;
	
	for( c = 0; c < cols; ++c ) {	  
	  if( c != r ) {	    
	    sumrow += trans[(r+(c*rows)) + k*(rows*cols)];	  
	  }	         
	}
      
	trans[(r+(r*rows)) + k*(rows*cols)] = (double)(1 - sumrow);       
      }        
    }
  }
}




