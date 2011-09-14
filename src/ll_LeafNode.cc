// $Id$
/**
 * @file
 * @brief The functions for the class LeafNode
 */
// $Log$

#include "ll.h"

// Bring in the global variables
extern int DIM;
extern double RADIUS;
extern double COMPACTNESS;

// -------------------------------------------------------- Constructor with parent provided
LeafNode::LeafNode (Node * parentNode ) : Node(parentNode )
{
  mdebug1("In LeafNode (Node * parentNode ) Constructor\n");
  m_children.push_back ( new Leaf );
  incrementNchilds();
}


// --------------------------------------------------------
bool LeafNode::copyTreeLeafNode( LeafNode * original, int kVar[], int kLeaf[], int *kLeafCount, int *leafCount, int *nrxL, int * oldDIM)
{
    mdebug4("\n In LeafNode() copyTreeLeafNode (START) \n");
    bool atLeast1Leaf = false;

    setNchilds(0);
    setNobs(0);
    while ( ((kLeaf[*kLeafCount]) <= (*leafCount) + original->getNchilds()) &  (*kLeafCount < *nrxL)){
        if (!atLeast1Leaf)  {
            m_children[0] = ( new Leaf(original->m_children[kLeaf[(*kLeafCount)] - (*leafCount) - 1], kVar ,oldDIM )); /*Copy leaf*/
            CF myCF = m_children[0]->getCF();
            setNobs(getNobs() + myCF.n);
            incrementNchilds();
            *kLeafCount = *kLeafCount + 1;
            atLeast1Leaf = true;
        }
        else {
            m_children.push_back ( new Leaf(original->m_children[kLeaf[(*kLeafCount)] - (*leafCount) - 1], kVar, oldDIM  )); /*Copy leaf*/
            CF myCF2 = m_children[m_children.size() - 1 ]->getCF();
            setNobs(getNobs() + myCF2.n);
            incrementNchilds();
            *kLeafCount = *kLeafCount + 1;
        }
    }
    *leafCount = (*leafCount) + original->getNchilds();
    recalculateCenter(); /*The selection of variables is handled by the Leaf copy contructor*/

    mdebug4("\n In LeafNode() copyTreeLeafNode (END) \n");
    return(atLeast1Leaf); //If we did not add leaf, we will delete this LeafNode after the return
}


// -------------------------------------------------------- Default constructor
LeafNode::LeafNode () : Node( NULL )
{
  mdebug1("In LeafNode () Constructor\n");
  m_children.push_back ( new Leaf() );
  incrementNchilds();
}

// -------------------------------------------------------- Constructor with parent & child provided
LeafNode::LeafNode (Node * parentNode, Leaf * newLeaf) : Node( parentNode )
{
  mdebug1("In LeafNode (Node * parentNode, Leaf * newLeaf) Constructor\n");

  m_children.push_back ( newLeaf );
  setNchilds(1);
  CF myCF = (*newLeaf).getCF();

  #ifdef debug4
    std::cout << "mydataDIM " << DIM;
  #endif

  //double newCenter[DIM]; variably-dimensioned array, specific to GNU
  double * newCenter = new double[DIM];
  
  for (int i=0; i < DIM; i++)
    newCenter[i] = myCF.sumXi[i] / myCF.n;
  updateCenter(newCenter, myCF.n);
}


// --------------------------------------------------------
LeafNode::~LeafNode() {
  mdebug1("In LeafNode()  Destructor\n");
  while (getNchilds() > 0){
    delete m_children[ getNchilds()-1 ];
    decrementNchilds();
  }
}




// --------------------------------------------------------
void LeafNode::addToNode(double data[], int obsnumber){
  mdebug1("In LeafNode::addToNode(double data[], int obsnumber)\n");

  // Find the closest Leaf
  int closest;   double newradius;
  findClosest ( data , closest, newradius);

  // Check condition
  CF newCF = m_children[closest]->CalcNewCF( data );
  double newcompact = CalcNewRadius ( newCF );

  if ( newcompact < COMPACTNESS &&  newradius < RADIUS ){
    m_children[closest]->addToLeaf ( newCF, newradius, obsnumber );
  }
  else {
    // Add a new leaf
    m_children.push_back (new Leaf( data, obsnumber));
    incrementNchilds();
  }
  updateCenter(data);

  if (getNchilds() > m_maxChildren)
    SplitNode();
}


// --------------------------------------------------------
double LeafNode::CalcNewRadius( CF & newData ){                                                   mdebug1("In LeafNode::CalcNewRadius( CF & newData )\n");

  double radius=0;
  for (int i=0; i< DIM; i++ )
    radius += newData.sumXisq[i + DIM*i]/ newData.n - square ( newData.sumXi[i] / newData.n );
  return ( radius );
}


// --------------------------------------------------------
void LeafNode::getLeaves(vector<Leaf*> & returndata){                                            mdebug1("In LeafNode::getLeaves\n");

  for (int i=0; i < getNchilds(); i++)
    returndata.push_back( m_children[i] );
}

// --------------------------------------------------------
void LeafNode::SplitNode(){                                                                      mdebug1("In LeafNode::splitNode\n");
                                                                                                 mdebug2("Splitting Leaf Node\n");

  unsigned int i;

    // Returns index of ones to go into new node
    vector<int> nodeGo = calcSplit();

    // Split this node into itself, and the new Node
    // Do it for the first one (then loop the rest)
    LeafNode * tmp = new LeafNode(getParent(), m_children[ nodeGo[0] ]);
    if (nodeGo.size() > 1){
      for (i=1; i < nodeGo.size(); i++){
    tmp->m_children.push_back(m_children[ nodeGo[i] ]);
    tmp->incrementNchilds();
      }
    }
    // Then remove from existing node
    for (i = nodeGo.size(); i > 0; i--){
      m_children.erase(m_children.begin() + nodeGo[i-1]);
      decrementNchilds();
    }

    // Recalculate each center
    tmp->recalculateCenter(); // New node
    recalculateCenter(); // Existing node

    // Give it to parent, and head on back up
    getParent()->addToNode(tmp);
}


// --------------------------------------------------------
double LeafNode::calcDistances(int i, int j){                                                    mdebug1("In LeafNode::calcDistances(int i, int j)\n");

  // Calculates the distance between leaves i and j
  double distance=0;

  CF CF1 = m_children[i]->getCF();
  CF CF2 = m_children[j]->getCF();

  for (int k=0; k < DIM; k++ )
    distance += square ( CF1.sumXi[k]/CF1.n - CF2.sumXi[k]/CF2.n );
  return ( distance );
}


// --------------------------------------------------------
void LeafNode::recalculateCenter(){                                                               mdebug1("In LeafNode::recalculateDistances\n");

  double newCenter[MAXD];
  int newN=0;

  // Initialise newCenter;
  for (int j=0; j < DIM; j++) newCenter[j] = 0;

  for (int i=0; i < getNchilds(); i++){
    CF myCF = m_children[i]->getCF();
    for (int j=0; j < DIM; j++)
      newCenter[j] += myCF.sumXi[j];
    newN += myCF.n;
  }

  for (int i=0; i < DIM; i++)
    newCenter[i] /= newN;
  updateCenter(newCenter, newN);
}


// --------------------------------------------------------
void LeafNode::findClosest( double data[], int & closest, double & rad){                            mdebug1("In Node::findClosest( double data[] )\n");

  int smallest=0;
  double tmpDistance, smDistance = m_children[0]->CalcDistance(data);

  // find the closest child node to this data point
  for (int i=1; i < getNchilds(); i++){
    tmpDistance = m_children[i]->CalcDistance(data);
    if (tmpDistance < smDistance){
      smallest = i;
      smDistance = tmpDistance;
    }
  }
  closest = smallest;
  rad = smDistance;
}
