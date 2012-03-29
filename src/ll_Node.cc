// $Id$
/**
 * @file
 * @brief The functions for the class Node
 */
// $Log$

#include "ll.h"
#include <typeinfo>

// Bring in the global variables
extern int DIM;
extern double RADIUS;

// --------------------------------------------------------
Node::Node()
{
    // Constructor: This means we are the parent
    mdebug1("In Node() Constructor\n");

    m_parentNode = NULL;
    m_nobs = 0;
    for (int i =0; i < DIM; i++) m_center[i] = 0;

    // And create a Leaf Node
    m_children.push_back(new LeafNode(this));
    m_nchilds = 1;

}


// --------------------------------------------------------
Node::Node( Node * parentNode )
{
    // Constructor: This means another node has created us
    mdebug1("In Node(Node * parentNode ) Constructor\n");

    m_nobs = 0;
    m_nchilds = 0;
    for (int i =0; i < DIM; i++) m_center[i] = 0;
    m_parentNode = parentNode;
}


// --------------------------------------------------------
bool Node::copyTree( Node * original, int kVar[], int kLeaf[], int *kLeafCount, int *leafCount, int *nrxL, int * oldDIM )
{
    mdebug4("\n \n In Node::copyTree(...) ");

    bool atLeast1Child = false;
    bool toAdd = false;  //Indicate if the new kid have at least 1 kid

    m_nchilds = 0;
    m_nobs = 0;

    for (int i = 0; i < original->getNchilds(); i++ )
    {
        if (original->m_children[i]->isLeafNode())  /* Children are LeafNodes*/
        {
            LeafNode * tmp = new LeafNode(this);
            toAdd = tmp->copyTreeLeafNode(dynamic_cast<LeafNode*> (original->m_children[i]),kVar , kLeaf, kLeafCount, leafCount , nrxL, oldDIM);
            if (toAdd == true)
            {
                m_children.push_back(dynamic_cast <Node*> (tmp));
                m_nchilds  =  m_nchilds +1;
                m_nobs = m_nobs + tmp->m_nobs;
            }
        }
        else   /* Children are Nodes*/
        {
            Node * tmp = new Node(this);
            toAdd = tmp->copyTree(dynamic_cast<Node*> (original->m_children[i]),kVar , kLeaf, kLeafCount, leafCount,  nrxL, oldDIM);
            if (toAdd == true)
            {
                m_children.push_back(tmp);
                m_nchilds  =  m_nchilds +1;
                m_nobs = m_nobs + tmp->m_nobs;
            }
        }
    } //End for
    recalculateCenter();

    if (m_nchilds != 0) atLeast1Child = true;
    return(atLeast1Child); //If new Node does not have at least one kid, we delete it after the return
}

// --------------------------------------------------------
Node::Node (Node * parentNode, Node * newNode)
{
    mdebug1("In Node (Node * parentNode, Node * newNode) Constructor\n");

    newNode->m_parentNode = this;
    m_parentNode = parentNode;
    m_children.push_back ( newNode );
    m_nchilds=1;
    updateCenter(newNode->m_center, newNode->m_nobs);
}

// --------------------------------------------------------

Node::~Node()
{
    // Destuctor
    mdebug1("In Node()  Destructor\n");

    while (m_nchilds > 0)
    {
        delete m_children[ m_nchilds-1 ];
        m_nchilds--;
    }
}


// --------------------------------------------------------
void Node::addToNode(double data[], int obsnumber)
{
    mdebug1("In Node::addToNode(double data[], int obsnumber)\n");
    mdebug1("In Node::addToNode(double data[], int obsnumber)\n" );

    // Find the nearest child node
    int nearest = findClosest(data);

    // Add to this node
    m_children[nearest]->addToNode( data , obsnumber);

    // Then update the current center stats for this node
    updateCenter(data);
}

// --------------------------------------------------------
void Node::addToNode(Node * newNode)
{
    mdebug1("In Node::addToNode (Node * newNode) \n");

    // Add new child
    m_children.push_back(newNode);
    m_nchilds++;

    // Do I need to split?
    if (m_nchilds > m_maxChildren)
        SplitNode();
}


// --------------------------------------------------------
double Node::CalcDistance(double data[])
{
    mdebug1("In Node::CalcDistance(double data[])\n");
    // Calculate the distance from data to my center

    double radius = 0;
    for (int i=0; i < DIM; i++) radius += square(data[i] - m_center[i]);
    return(radius);
}

// --------------------------------------------------------
int Node::findClosest( double data[] )
{
    mdebug1("In Node::findClosest( double data[] )\n");

    int smallest=0;
    double tmpDistance, smDistance = m_children[0]->CalcDistance(data);
#ifdef debug3
//    std::cout << "Distances 1:\n" << smDistance;
#endif

    // find the closest child node to this data point
    for (int i=1; i < getNchilds(); i++)
    {
        tmpDistance = m_children[i]->CalcDistance(data);
#ifdef debug3
//        std::cout << " " << tmpDistance;
#endif
        if (tmpDistance < smDistance)
        {
            smallest = i;
            smDistance = tmpDistance;
        }
    }
    return( smallest );
}

// --------------------------------------------------------
void Node::returnData(SEXP * output)
{
    mdebug1("In Node::returnData\n");

    vector <Leaf *> returndata;

    getLeaves(returndata);
    int nleafs = returndata.size();

    SEXP NinLeaf, LeafMembers, SumXi, SumXiSq;
    PROTECT ( NinLeaf = allocVector ( INTSXP, nleafs ) );
    PROTECT ( LeafMembers = allocVector ( VECSXP, nleafs ) );
    PROTECT ( SumXi = allocMatrix ( REALSXP, DIM, nleafs ) );
    PROTECT ( SumXiSq = alloc3DArray ( REALSXP, DIM, DIM, nleafs ) );

    for (int i=0; i < nleafs; i++)
    {
        returndata[i]->returnData ( &NinLeaf, &LeafMembers,
                                    &SumXi, &SumXiSq, i);
    }
    SET_VECTOR_ELT ( *output, 0, NinLeaf );     /* N */
    SET_VECTOR_ELT ( *output, 1, SumXi );       /* sumXi */
    SET_VECTOR_ELT ( *output, 2, SumXiSq );     /* sumXisq */
    SET_VECTOR_ELT ( *output, 3, LeafMembers ); /* members */
    UNPROTECT ( 4 );
}

// --------------------------------------------------------
void Node::getLeaves(vector<Leaf *> & returndata)
{
    mdebug1("In Node::getLeaves\n");

    for (int i=0; i < m_nchilds; i++)
        m_children[i]->getLeaves(returndata);
}

// --------------------------------------------------------
void Node::updateCenter(double data[])
{
    mdebug1("In Node::updateCenter(double data[])\n");

    for (int i=0; i<DIM; i++)
        m_center[i] = (m_center[i]*m_nobs + data[i])/(m_nobs+1);
    m_nobs++;
}

// --------------------------------------------------------
void Node::updateCenter(double data[], int nobs)
{
    mdebug1("In Node::updateCenter(double data[], int nobs)\n");

    for (int i=0; i<DIM; i++)
        m_center[i] = data[i];
    m_nobs = nobs;
}



// --------------------------------------------------------
void Node::SplitNode()
{
    mdebug1("In Node::SplitNode\n");

    // Check first if parent
    if (m_parentNode == NULL)
    {

        mdebug2("Splitting Parent Node\n");

        // OK, I'm the parent - create a new node
        Node * tmp = new Node(this);

        // Tell the children who's the new parent
        for (int i=0; i < m_nchilds; i++)
            m_children[i]->m_parentNode = tmp;

        tmp->updateCenter(m_center, m_nobs);
        tmp->m_children = m_children;
        tmp->m_nchilds = m_nchilds;

        // Now wipe my slate clean
        m_children.clear();
        m_children.push_back(tmp);
        m_nchilds = 1;
    }
    else
    {
        mdebug2("Splitting Node\n");

        // Returns index of ones to go into new node
        vector<int> nodeGo = calcSplit();

        // Split this node into itself, and the new Node
        // Do it for the first one (then loop the rest)
        Node * tmp = new Node(m_parentNode, m_children[ nodeGo[0] ]);
        unsigned int i;
        if (nodeGo.size() > 1)
        {
            for (i=1; i < nodeGo.size(); i++)
            {
                tmp->m_children.push_back(m_children[ nodeGo[i] ]);
                tmp->m_nchilds++;
            }
        }

        // Then remove from existing node
        for (i = nodeGo.size(); i > 0; i--)
        {
            m_children.erase(m_children.begin() + nodeGo[i-1]);
            m_nchilds--;
        }

        // update center
        tmp->recalculateCenter(); // New node
        recalculateCenter(); // Existing node

        // Give it to parent, and head on back up
        getParent()->addToNode(tmp);
    }
}

// --------------------------------------------------------
vector<int> Node::calcSplit()
{
    // nodeA is the one thats being removed from the exisiting leaf
    mdebug1("In Node::calcSplit\n");

    int i,j;
    //double left, right; //
    int nchilds = getNchilds();
    vector<int> nodeA;

    //double distances[nchilds*nchilds]; variably-dimensioned array, specific to GNU
    double * distances = new double[nchilds*nchilds];

    //for (i=0; i < nchilds; i++)
    //    for (j=0; j < nchilds; j++)
    //        distances[i + nchilds*j] = calcDistances(i,j); //LYS
    //        distances[i + nchilds*j] = calcDistances(i,j); //LYS
    distances[0] = 0;
    for ( i=0; i < nchilds-1; i++)
        for ( j=i+1; j < nchilds; j++){
            distances[i + nchilds*j] = calcDistances(i,j);  }
    // NB Could be more efficient here (with a upper-diag matrix) OK DONE

    vector<int> widest = which_min(distances, nchilds, FALSE);

    for (i=0; i < widest[0]; i++)
    {
      if( distances[ i + nchilds*widest[1] ] < distances[ i + nchilds*widest[0] ] )
        nodeA.push_back(i);
    }
     for (i=widest[0]; i < widest[1]; i++)
    {
         if( distances[ i + nchilds*widest[1] ] < distances[ widest[0] + nchilds*i ] )
            nodeA.push_back(i);
    }
     for (i=widest[1]; i < nchilds-1; i++)
    {
        if( distances[ widest[1] + nchilds*i ] < distances[ widest[0] + nchilds*i ] )
            nodeA.push_back(i);
    }


    return (nodeA);
}

// --------------------------------------------------------
double Node::calcDistances(int i, int j)
{
    // Calculates the distance between leaves i and j
    mdebug1("In Node::calcDistances(int i, int j)\n");

    double distance=0;
    for (int k=0; k < DIM; k++ )
        distance += square ( m_children[i]->m_center[k] -  m_children[j]->m_center[k]);
    return ( distance );
}

// --------------------------------------------------------
void Node::recalculateCenter()
{
    mdebug1("In Node::recalculateCenter\n");
    double newCenter[MAXD];
    int newN=0;

    // Initialise newCenter;
    for (int j=0; j < DIM; j++) newCenter[j] = 0;
    for (int i=0; i < m_nchilds; i++)
    {
        // Get pointer
        Node * child = m_children[i];
        for (int j=0; j < DIM; j++) newCenter[j] += child->m_center[j] * child->m_nobs;
        newN += child->m_nobs;
    }
    for (int i=0; i < DIM; i++)
        newCenter[i] /= newN;
    updateCenter(newCenter, newN);
}
