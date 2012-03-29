#ifndef LL_H_
#define LL_H_
#define MAXD 30 /* Specify max number of dimensions */
#define MAXL 40 /* Max number of Node Children */
#define MAXB 40 /* Max number of LeafNodeChildren */

//#define debug /* debug flag */
//#define debug2 /* Minimal */
//#define debug3 /* for leaf decisions */
//#define debug4  /*For variably-dimensioned arrays and birch.[*/

#include <iostream>
#include <Rinternals.h> /* for REAL, SEXP etc */
#include <vector>
#include <cmath>
#include <typeinfo>
using std::vector;

#define CHECK_BIRCH_OBJECT(s) do { \
    if (TYPEOF(s) != EXTPTRSXP || \
        R_ExternalPtrTag(s) !=  RRP_birch_tag || \
    R_ExternalPtrAddr ( XPtr ) == NULL) \
        error("bad birch object"); \
} while (0)


//-------------------------- Leaf Class
//! The CF Struct
/*! This struct contains the actual clustering features (summary statistics) */
  struct CF {
    int n;            /*!<The number of observations */
    double sumXi[MAXD];           /*!<The cumsum of the observations */
    double sumXisq[MAXD * MAXD];  /*!< The cumsum squared observations */
  };

//! The Leaf Class
/*!
  Contains the actual leaves (subclusters), which are eventually returned.
*/
class Leaf {
 private:
  CF m_CF; /*!< The Clustering Feature struct */
  double m_radius;                        /*!< The radius of the leaf */
  vector<int> m_observations;             /*!< Vector of memberships */

 public:
  // Constructors and Destructors
  //! The main constructor (with data attached)
  Leaf (double * data, int obsnumber);
  //! Constructor for copy (with original provided)
  Leaf ( Leaf * original, int kVar[], int * oldDIM );

  //! The default constructor
  Leaf ();

  //! The default destructor
  ~Leaf ();

  // Other functions
  //! Adds data to the leaf
  void addToLeaf ( CF newCF, /*!<  updated CF with new data */
           double newRadius, /*!< the new Radius after data added */
           int obsnumber /*!< which observation number */
           );

  //! Calculates distance to center
  double CalcDistance ( double * y /*!< observation from which to calc dist to center */
            );

  //! Calculates a new CF, based on data given to it
  CF CalcNewCF ( const double * data /*!< data to add to CF */
         );

  //! return the data from the leaf
  /*! Passed a bunch of SEXP types to populate */
  void returnData ( SEXP * NinLeaf, SEXP * LeafMembers,
            SEXP * SumXi, SEXP * SumXiSq, int index);

  //! Returns the CF of the leaf
  CF getCF() { return m_CF;}

  //! Returns the radius of the leaf
  double getRadius() { return m_radius; }

  //! Returns the number of observations in the leaf
  int getN() {return m_CF.n; }
};



//-------------------------- Node Class
//! The Node Class
/*!
  This is one of the most important classes - data is passed to it (the parent Node), who then decides how it is passed down the tree.
  The main go to guy!
*/

class Node {
 private:
  Node * m_parentNode; /*<! A pointer to the parent (if not the parent) */
  static const int m_maxChildren = MAXL; /*<! how many children before splitting*/
  vector<Node *> m_children; /*<! a STL vector containing pointers to my children */
  int m_nchilds; /*<! how many children do I have */
  double m_center[MAXD]; /*<! my center */
  int m_nobs; /*<! how many observations in total below me */

  //! Let's split myself
  void SplitNode();

  //! Calcs the distance between my centerand data
  double CalcDistance(double * data /*!< object to calc distance to */
              );
  //! Recalculate my center
  void recalculateCenter();

  //! Find the closest child to pass down to
  int findClosest(double * data /*<! data looking for a home */
          );

 protected:
  //! Passing a single data point, and update my center
  void updateCenter(double * data /*<! the data being passed */
            );

  //! Passing a bunch of data points, and update my center
  void updateCenter(double * data, /*<! the data being passed */
            int nobs /*<! the number of observations being passed */
            );

  //! If splitting, tell me which children go where
  virtual vector<int> calcSplit();

  //! Calculate the distance between two children
  virtual double calcDistances(int i, int j);

  // Accessor functions for LeafNode
  //! incrementNchilds
  void incrementNchilds() {m_nchilds++;}

  //!decrementNchilds
  void decrementNchilds() {m_nchilds--;}

  //!Specify Nchilds
  void setNchilds(int nchilds){m_nchilds = nchilds;}

  //!Specify Nobs
  void setNobs(int nobs){m_nobs = nobs;}

  //! Who's my daddy?
  Node * getParent() {return m_parentNode;}


 public:

  //! Default constructor
  Node();
  //! Default constructor without empty LeafNode
  //Node(int ind, int ind2);
  //! Constructor with parent information
  Node( Node * parentNode );
  //! Constructor with parent and child information
  Node (Node * parentNode, Node * newNode);

  bool copyTree( Node * original, int kVar[], int kLeaf[], int *kLeafCount, int *leafCount, int *nrxL, int * oldDIM );

  //To be use with addChildren
  //void addChildren(Node * tmp) {m_children.push_back( tmp) ;}
  //To use incrementNchilds wich is protected
  //void increment(){incrementNchilds() ;}
  //Node getChild(int i) {return(m_children[i]) ;}

  virtual bool isLeafNode(){ /*std::cout << "\n Node::IsLeafNode " ;*/ return(false) ;}

  //! Return Nchilds
  int getNchilds() {return m_nchilds;}

  //!Return Nobs
  int getNobs(){return m_nobs ;}

  //! Default destructor
  virtual ~Node();
  //! Add data to this node
  /*!
     The workhorse of the entire algorithm. Decides who's best, then passes it down (often to itself in a Node class lower).
   */
  virtual void addToNode(double * data, int obsnumber);
  //! A (sub) node to my node. For splitting
  virtual void addToNode(Node * newNode);
  //! Part of the process of getting a STL vector pointers to leaves (for returning the data)
  virtual void getLeaves(vector<Leaf *> & returndata);
  //! Returning the data
  void returnData(SEXP * output);
};



//-------------------------- LeafNode Class
//! The LeafNode class
/*!
  The class that interfaces between Nodes above and Leafs below. Inherits from Node.
*/
class LeafNode : public Node // Inherits from Node class
{
 private:
  static const int m_maxChildren = MAXB; /*<! How many children am I allowed? */
  vector<Leaf *> m_children; /*<! A STL vector of pointers to my children */

  //! Given a CF, what would my new radius be?
  double CalcNewRadius( CF & newData );
  //! Split myself
  void SplitNode();
  //! recalculate my center based on information already there
  void recalculateCenter();
  //! find my closest child
  void findClosest(double * data, int & closest, double & rad);

 protected:
  //! Part of the process of returning a vector of pointers to the leafs
  virtual void getLeaves(vector<Leaf *> & returndata);
  //! calculate the distance between leaf i and j (part of the splitting process)
  virtual double calcDistances(int i, int j);

 public:

  //! Constructor (with parent provided)
  LeafNode(Node * parentNode);
  //! Constructor (with parent and child provided)
  LeafNode (Node * parentNode, Leaf * newLeaf);

  bool copyTreeLeafNode( LeafNode * original, int kVar[], int kLeaf[], int *kLeafCount, int *leafCount, int *nrxL, int * oldDIM );
  virtual bool isLeafNode(){/*std::cout << "\n LeafNode::IsLeafNode " ;*/ return(true) ;}

  //! Default constructor
  LeafNode();
  //! Default destructor
  ~LeafNode();
  //! add data (with obsnumber) to leaf below
  virtual void addToNode(double * data, int obsnumber);
};



// Inline functions --------------------------------------
inline double square ( double x ) { return x * x; }


inline void mdebug1 ( const char output[] ) {
#ifdef debug
//  std::cout << output;
#endif
}
inline void mdebug2 ( const char output[] ) {
#ifdef debug2
  //std::cout << output;
#endif
}
inline void mdebug3 ( const char output[] ) {
#ifdef debug3
  //std::cout << output;
#endif
}
inline void mdebug4 ( const char output[] ) {
#ifdef debug4
 // std::cout << output;
#endif
}



// Other definitions
vector<int> which_min(double distances[], int nrx, bool min=TRUE);

#endif
