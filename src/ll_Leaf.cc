// $Id$
/**
 * @file
 * @brief The functions for the class Leaf
 */
// $Log$

#include "ll.h"
// Bring in the global variables
extern int DIM;
extern double RADIUS;

// --------------------------------------------------------

// constructors
Leaf::Leaf ( double data[], int obsnumber ) // Constructor with data
{
  mdebug1("In Leaf( double data[], int obsnumber ) Constructor\n");

    int i,j;
    for ( i=0; i < DIM; i++ ) m_CF.sumXi[i] = data[i];

    for ( i=0; i < DIM; i++ )
        for ( j=0; j < DIM; j++ )
            m_CF.sumXisq[i + DIM*j] = data[i] * data[j];

    m_radius = 0;
    m_CF.n = 1;
    m_observations.push_back ( obsnumber );
}


// --------------------------------------------------------
Leaf::Leaf () { // Constructor without data

  mdebug1("In Leaf()  Constructor\n");

// Initialize CF
    for ( int i=0; i < DIM; i++ ) {
        m_CF.sumXi[i] = 0;
        for ( int j=0; j < DIM; j++ )
            m_CF.sumXisq[i + DIM*j] = 0;
    }
// Initialize radius and n
    m_CF.n = 0;
    m_radius = 0;
}

 

// --------------------------------------------------------
// destructors
Leaf::~Leaf() { mdebug1("In Leaf()  Destructor\n"); }

// other methods
// --------------------------------------------------------
void Leaf::addToLeaf ( CF newCF, double newRadius, int obsnumber ) {

  mdebug1("In Leaf::addToLeaf \n");

    // Add new observation number
    m_observations.push_back ( obsnumber );

    // Update CF
    m_CF = newCF;

    // Calculate the new radius
    m_radius = newRadius;
}	



// --------------------------------------------------------
double Leaf::CalcDistance ( double y[] ) {
  mdebug1("In Leaf::CalcDistance\n");
 
  double distance=0;
  if (m_CF.n != 0) {
    for (int i=0; i < DIM; i++ )
      distance += square ( m_CF.sumXi[i]/m_CF.n - y[i] );
  }
  return ( distance );
}



// --------------------------------------------------------
CF Leaf::CalcNewCF ( const double data[] ){
  mdebug1("In Leaf::CalcNewCF\n");

 int i, j;

 CF newData = m_CF;

  /* Calculate new sumXi and sumXisq */
  for ( i=0; i < DIM; i++ ) newData.sumXi[i] +=  + data[i];

  for ( i=0; i < DIM; i++ )
    for ( j=0; j < DIM; j++ )
      newData.sumXisq[i + DIM*j] += data[i] * data[j];

  newData.n++;
  return(newData);
}


// --------------------------------------------------------
void Leaf::returnData ( SEXP * NinLeaf, SEXP * LeafMembers, 
			SEXP * SumXi, SEXP * SumXiSq, int index)
{
  int j, *INinLeaf, *ILeafMemTmp;
  double *RSumXi, *RSumXiSq;

  // Number in leaf
  INinLeaf = INTEGER ( *NinLeaf );
  INinLeaf[index] = m_CF.n;
  
  // Leaf members
  SEXP LeafMemTmp;
  PROTECT ( LeafMemTmp = allocVector ( INTSXP, m_CF.n ) );
  ILeafMemTmp = INTEGER ( LeafMemTmp );
  for ( j=0; j < m_CF.n; j++ )
    ILeafMemTmp[j] = m_observations[j];
  SET_VECTOR_ELT ( *LeafMembers, index, LeafMemTmp ); 
  UNPROTECT(1);

  // sumXi
  RSumXi = REAL ( *SumXi );
  for ( j=0; j < DIM; j++ )
    RSumXi[index*DIM + j] = m_CF.sumXi[j];
  
  // sumXisq
  RSumXiSq = REAL ( *SumXiSq );
  for ( j=0; j < square ( DIM ); j++ )
    RSumXiSq[index*DIM*DIM + j] = m_CF.sumXisq[j];
      
}

