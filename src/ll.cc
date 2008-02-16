// $Id$
/**
 * @file
 * @brief The main API for R (suitably bound in C externs)
 * @author Me
 * @version 1.1-2
 */
// $Log$

#include "ll.h"
// Define global variables
int DIM;
double RADIUS;
double COMPACTNESS;
static SEXP RRP_birch_tag;
static void _birch_free(SEXP XPtr);

/**
   The function LL_main is the work horse. It initilizes the tree, takes in a matrix feeds the underlying object 
   headed up by Parent (a pointer to the class Node), and then returns the data, or the tree.
   
   @param[in] X the matrix containing data
   @param[in] R the radius criteria for forming the tree
   @param[in] COMPACT the compactness criteria for forming the tree
   @param[in] RETTREE either 1 or 0. If 1, then keep the tree in memory.
   @return outputs (a list, containing the usual CFs), pointer (of type SEXP, but this time contains a pointer instead (if RETTREE is 1))
*/
extern "C" {

  SEXP LL_main ( SEXP X, SEXP R, SEXP COMPACT, SEXP RETTREE) {
    // assign global variables 
    DIM = ncols( X );
    RADIUS = REAL(R)[0];
    COMPACTNESS = REAL(COMPACT)[0];
    
    //  Create the top node
    mdebug2("Creating Parent Node\n");

    Node * Parent = new Node;

    // Load the data
    int nrx = nrows ( X );
    double mydata[DIM];
    for (int i=0; i < nrx; i++ ) {
      for (int j=0; j < DIM; j++ ) 
	mydata[j] = REAL ( X ) [i + j*nrx];

#ifdef debug2
      std::cout << "\nAdding Observation " << i << std::endl;;
#endif
    
      Parent->addToNode ( mydata, i + 1 );
    }

    mdebug2("\nReturning Data\n");
    if (INTEGER(RETTREE)[0] == 1) {
      // Return a pointer
      SEXP pointer;
      PROTECT(pointer = R_MakeExternalPtr (Parent, RRP_birch_tag, R_NilValue));
      R_RegisterCFinalizerEx(pointer, _birch_free, TRUE);
      UNPROTECT(1);
      return ( pointer );
    }
    else {
      //  Return a list of memberships
      SEXP outputs;
      PROTECT ( outputs = allocVector ( VECSXP,  4) );
      Parent->returnData ( &outputs);
      delete Parent; 
      UNPROTECT ( 1 );
      return ( outputs );
    }
  }

} 

/**
   The function LL_killtree removes the tree from memory (by firing off a series of destructors), and then nulls the pointer. Returns NULL.
   @param[in] XPtr the external pointer from which to start the killing spree
*/
extern "C" {
    SEXP LL_killtree ( SEXP XPtr){
    CHECK_BIRCH_OBJECT ( XPtr );
    Node * Parent = (Node *) R_ExternalPtrAddr ( XPtr ) ;
    delete Parent;
    R_ClearExternalPtr(XPtr);
    return R_NilValue;
  }
}

void _birch_free(SEXP XPtr){
  Node * Parent = (Node *) R_ExternalPtrAddr ( XPtr ) ;
  delete Parent;
  R_ClearExternalPtr(XPtr);
}

/**
   The function LL_getdata recieves a pointer (which is an external pointer telling us where the parent is). 
   We then call returndata (part of the class Node), passing by assignment the SEXP type output vector
   @param[in] XPtr the external pointer from which to start the loading
   @return outputs (of type SEXP. Is a list, containing the usual CFs.)
*/
extern "C" {
  SEXP LL_getdata ( SEXP XPtr){
    CHECK_BIRCH_OBJECT ( XPtr );
    SEXP outputs;
    Node * Parent = (Node *) R_ExternalPtrAddr ( XPtr ) ;
    PROTECT ( outputs = allocVector ( VECSXP,  4) );
    Parent->returnData ( &outputs);
    UNPROTECT ( 1 );
    return ( outputs );
  }
}

/** 
    The function LL_adddata takes a pointer to the Parent, and a matrix SEXP, and loads more data onto the tree.
    It keeps the same radii as before, and returns NULL.
    @param[in] XPtr the external pointer from which to start the loading
    @param[in] X the matrix containing data
*/
extern "C" {
   SEXP LL_adddata ( SEXP XPtr , SEXP X) {
     CHECK_BIRCH_OBJECT ( XPtr );
     Node * Parent = (Node *) R_ExternalPtrAddr ( XPtr ) ;
     int nrx = nrows ( X );
     double mydata[DIM];
     for (int i=0; i < nrx; i++ ) {
       for (int j=0; j < DIM; j++ ) 
	 mydata[j] = REAL ( X ) [i + j*nrx];
       
#ifdef debug2
       std::cout << "\nAdding Observation " << i << std::endl;;
#endif
       Parent->addToNode ( mydata, i + 1 );
     }
     
     return R_NilValue;
   }
}

/**
   The function LL_getdim simply interrogates the tree for its dimensions (number of leaves, number of observations)
   Returns a vector containing these.
   @param[in] XPtr the external pointer from which to start the loading
   @return  outputs (SEXP vector of length two.)
*/
extern "C" {
  SEXP LL_getdim ( SEXP XPtr){
    CHECK_BIRCH_OBJECT ( XPtr );
    SEXP outputs;
    Node * Parent = (Node *) R_ExternalPtrAddr ( XPtr ) ;

    vector <Leaf *> returndata;
    Parent->getLeaves(returndata);
    int nleafs = returndata.size();
    
    int N =0;
    for (int i=0; i < nleafs; i++)
      N+= returndata[i]->getN();
    
    PROTECT ( outputs = allocVector (INTSXP ,  2) );
    INTEGER(outputs)[0] = nleafs;
    INTEGER(outputs)[1] = N;
    UNPROTECT ( 1 );
    return ( outputs );
  }
}
