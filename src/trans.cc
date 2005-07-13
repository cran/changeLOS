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
// .C("trans", as.character(state.names), 
//             as.integer(length(state.names)),
//             as.character(state.names[model$jumps[,1]]), 
//             as.character(state.names[model$jumps[,2]]), 
//             as.integer(nrow(model$jumps)),
//             nj = as.integer(nrjumps[,3]), 
//             as.integer(nr.start) , 
//             as.double(times), 
//             as.integer(length(times)),
//             as.character(observ$from), 
//             as.character(observ$to), 
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

#include "matrix.h"

#include <map>
#include <string>

using namespace std;



typedef map<string,int> StringIntMap;  // map with string keys and int values 
typedef StringIntMap::iterator itStringIntMap;

typedef map<double,int> DoubleIntMap;
typedef DoubleIntMap::iterator itDoubleIntMap; // map with double keys and int values

extern "C" {

  void trans( char** state_names, 
              int* nr_states, 
              char** from, 
              char** to, 
              int* nr_from_to, 
              int* nr_jumps, 
              int* nr_start, 
              double* times, 
              int* nr_times,
	      char** observ_from, 
              char** observ_to, 
              double* observ_time, 
              int* nr_observ,
              double* a, 
              int* nrow, 
              int* ncol, 
              int* n ) {

    StringIntMap StateNames;  // for holding the state names
    
    StringIntMap StatesFrom;  // for holding the states from where jumps are possible

    string  CensState = "";    // name of the censoring state

    DoubleIntMap Times;       // for holding the timepoints at which the jumps occurs

 

    // storing the passed state names
    for( int i = 0; i < *nr_states; ++i ) {
      StateNames.insert( StringIntMap::value_type(string(state_names[i]),i) );     
    }
    
    // the last passed statename is the name of the censoring variable   
    if( StateNames.size() > 0 ) {
      itStringIntMap pos = StateNames.end(); // position after the last element
      --pos;
      CensState = pos->first;       
    }     
    
    // storing the passed states from where jumps are possible
    for( int i = 0; i < *nr_from_to; ++i ) {
      nr_jumps[i] = 0;      
      if( StatesFrom.find( string(from[i]) ) == StatesFrom.end() ) {
	StatesFrom.insert( StringIntMap::value_type(string(from[i]),i) );        
      }
    }      
    
    // storing the passed timepoints at which jumps occurs       
    for( int i = 0; i < *nr_times; ++i ) {
      Times.insert( DoubleIntMap::value_type(times[i],i) );      
    }

    // initialize the array of the transition matrices
    Array A(a,*nrow,*ncol,*n);

    //cout << A << endl;

    // compute the number of jumps and store it in the array of the transition matrices
    for(int i=0 ; i < *nr_observ; ++i) {      // loop over the observations (jumps)
      string o_from(observ_from[i]);   // state name from where the observed jump occurs   
      string o_to(observ_to[i]);       // state name to which the observed jump occurs

      // add the jump to the total number of jumps from state 'o_from' to state 'o_to'
      for( int j = 0; j < *nr_from_to; ++j ) {
	if( o_from == string(from[j]) && o_to == string(to[j]) ) {
	  nr_jumps[j] += 1;
	}
      }

      if( o_to != CensState ) {   // jumps to the censoring state (censoring in state 'o_from') are ignored
        // find out the observed time point
	itDoubleIntMap pos = Times.find( observ_time[i] );

        if( pos != Times.end() ) {
          // now we have the observed time point k
	  int k = pos->second;    

          // find out the state from where the observed jump occurs
	  itStringIntMap posfrom = StateNames.find(o_from);

          if( posfrom != StateNames.end() ) {
            // now we have the state (the position r in the vector of statenames) from where the observed jump ocurs
	    int r = posfrom->second;

            // find out the state to which the observed jump occurs
	    itStringIntMap posto = StateNames.find(o_to);
            
	    if( posto != StateNames.end() ) {
              // now we have the state (the position c in the vector of statenames) to which the observed jump ocurs
	      int c = posto->second;
              
              // add this jump in the k'th matrix at column c and row r
	      A[k][r][c] += 1;
	    }
	  }	  
	}	
      }
    }

    //cout << A << endl;

    // compute the risk
    Array B(A);
    
    for(int i=0; i < *nr_times; ++i) {     
      for( itStringIntMap frompos = StatesFrom.begin(); frompos != StatesFrom.end(); ++frompos ) { 
        itStringIntMap pos =  StateNames.find(frompos->first);;
        int j = pos->second;        
	double risk = nr_start[j];        // Number in state j at start        
	double jumps_to_j_before_i = 0;   // jumps to state j before timepoint number i
        double jumps_from_j_before_i = 0; // jumps from state j before timepoint number i
	double cens_from_j_before_i = 0;  // censoring in state j before timepoint number i

	for( int k = 0; k < i; ++k ) {
	   for( int r = 0; r < *nrow; ++r ) {
	     jumps_to_j_before_i += B[k][r][j];
	   }

	   for( int c = 0; c < *ncol; ++c ) {
	     jumps_from_j_before_i += B[k][j][c];
	   }
        }

        for(int ii = 0 ; ii < *nr_observ; ++ii) {
          string o_from(observ_from[ii]);
          string o_to(observ_to[ii]);

          if( o_to == CensState ) {                        
            if( observ_time[ii] < times[i] ) {
	      itStringIntMap posfrom_1 = StateNames.find(o_from);
	      
              if( posfrom_1 != StateNames.end() && posfrom_1->second == j ) {
	        cens_from_j_before_i += 1;
	      }	  
	    }	
          }
        }
	
	risk = risk + jumps_to_j_before_i - jumps_from_j_before_i - cens_from_j_before_i;
    
	if( risk > 0 ) {
	   for( int c = 0; c < *ncol; ++c ) {             
	     A[i][j][c] = A[i][j][c] / risk;
	   }
	}  
      }
     
      // compute the diagonal elements: the sum over each row must be 1
      for( int r = 0; r < *nrow; ++r ) {
        double sumrow = 0;

        for( int c = 0; c < *ncol; ++c ) {
	  if( c != r ) {
	    sumrow += A[i][r][c];
	  }	 
        }

        A[i][r][r] = (double)(1 - sumrow);
      }               
    }
    
    //cout << A << endl;

    // prepare the array of the transition matrices for returning (transform the array to an vector)
    A.as_double(a);
  }
}
