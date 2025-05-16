#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <Eigen/Core> // Core Eigen functionality
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <lattice.h>


template <typename T>
using Tensor2D = std::vector<std::vector<T>>;

class Node {
public:
    double data;
    int spin; 
    Node* next;

    Node(double val) {
        data = val;
        spin = 0; 
        next = nullptr;
    }
};


class AdjacencyMatrix {
private:
     std::vector<std::vector<bool>> adjacency;
public:
    AdjacencyMatrix(int n) : adjacency(n, std::vector<bool>(n, 0)) {}

    bool isConnected(int i, int j)
    {
        if (i >= 0 && i < adjacency.size() && j >= 0 && j < adjacency.size()) {
        return adjacency[i][j];}
        else{return false;}
    }
    void AddLink(int i, int j)
    {
        adjacency[i][j] = true; 
        adjacency[j][i] = true; 
    }

    void DeleteLink(int i, int j)
    {
        adjacency[i][j] = false; 
        adjacency[j][i] = false; 
    }
};
class UnionFind {
public:
    std::vector<int> parent;
    std::vector<int> rank;
    std::unordered_map<Node*, int> nodeToId;
    std::unordered_map<int, Node*> idToNode;
    
    int nextId = 0;
    // Constructor: Initialize the data structure
    UnionFind() {}; 

    // Add a Node to the Union-Find structure
    void addNode(Node* node) {
        if (nodeToId.find(node) == nodeToId.end()) {
            int id = nextId++;
            nodeToId[node] = id;
            parent.push_back(id);
            idToNode[id] = node;
            rank.push_back(0);  // Initial rank is 0 for new nodes
        }
    }

    // Find the root of the set containing the node
    int find(Node* node) {
        if (nodeToId.find(node) == nodeToId.end()) {
            throw std::runtime_error("Node not found in Union-Find structure");
        }
        int id = nodeToId[node];

        while (parent[id] != id) {
            parent[id] = parent[parent[id]];  // Path compression
            id = parent[id];
        }
        return id;
    }
    // Union the sets containing nodes x and y
    void unionSets(Node* x, Node* y) {
        int rootX = find(x);
        int rootY = find(y);

        if (rootX != rootY) {
            if (rank[rootX] > rank[rootY]) {
                parent[rootY] = rootX;
            } else if (rank[rootX] < rank[rootY]) {
                parent[rootX] = rootY;
            } else {
                parent[rootY] = rootX;
                rank[rootX] += 1;
            }
        }
    }

    // Check if two nodes are in the same set
    bool connected(Node* x, Node* y) {
        return find(x) == find(y);
    }

    std::vector<Node*> FindComponents(Node* node)
    {
        std::vector<Node*> res; 
        for (int i = 0; i<nextId; i++)
        {
            if (find(node)==find(idToNode[i]))
            {res.push_back(idToNode[i]); }    
        }
        return res; 
    }

    void PrintComponents()
    {
        for (int i = 0; i < nextId; i++) {
            Node* node = idToNode[i];
            int root = find(node);
            std::cout << "Node ID " << i << " is in component with root ID " << root << std::endl;
        }
    }

    void reset() {
        parent.clear();
        rank.clear();
        nodeToId.clear(); 
        idToNode.clear(); 
        nextId=0; 
    }
};


class CircularLinkedList {
public:
    double invT; 
    int member_num; 
    Node* head; 
    CircularLinkedList(double val) {
        invT = val; 
        member_num = 0; 
        head = nullptr;
    }


      ~CircularLinkedList() {
        Node* current = head;
        if (current) {
            Node* start = head;
            do {
                Node* next = current->next;
                current->next = nullptr;
                delete current;
                current = next;
            } while (current != start);
        }
    }

    void insertAfter(Node* prevNode, double value) {
        if (prevNode == nullptr) return;

        Node* newNode = new Node(value);
        newNode->next = prevNode->next;
        prevNode->next = newNode;
        member_num ++; 
    }

    void insertAtEnd(double value) {
        Node* newNode = new Node(value);
        if (head == nullptr) {
            head = newNode;
            head->next = head;
            
        } else {
            Node* temp = head;
            while (temp->next != head) {
                temp = temp->next;
            }
            temp->next = newNode;
            newNode->next = head;
        }
        member_num ++;
    }

    void placeSpin(Node* Node, int Spin){
        Node->spin = Spin; 
    }

    void deleteNode(Node* node){
        if (node == nullptr||head==nullptr) return;
        if (node == head) {
            // Find the node before head
            Node* temp = head;
            while (temp->next != head) {
                temp = temp->next;
            }
            // Only one node in the list
            if (temp == head) {
                delete head;
                head = nullptr;
            } else {
                Node* newHead = head->next;
                temp->next = newHead;
                delete head;
                head = newHead;
            }
        }

        else {
            Node* prev = head;
            while (prev->next != head && prev->next != node) {
                prev = prev->next;
            }
            // Node not found
            if (prev->next != node) return;

            prev->next = node->next;
            delete node;
        }
        member_num--;
    }
    double findInterval(Node* prevNode){//head 
        double val1  = prevNode->data;
        double val2  = prevNode->next->data;

        if (prevNode==head&&head->next == head){return invT;}
        else if (val2>val1){return val2-val1;}
        else{return invT+(val2-val1);} //Wrap around has occurred. 
    }
    std::vector<double> returnAllIntervals()
    {
        std::vector<double> intervals;
        Node* temp = head; 
        if (temp==nullptr){return intervals;}
        do {
        double dt = findInterval(temp);
        intervals.push_back(dt);
        temp = temp->next;
        } while (temp != head);
        return intervals;
    }
    void mergeDuplicateSpins()
    {
        if (!head|| head->next==head){return;}
        Node* temp = head;
        do{     
        while(temp->next->spin == temp->spin)
        {deleteNode(temp->next);};
        temp = temp->next;
        }while(temp!=head);
    }

    void display() {
    if (head == nullptr) return;

    Node* temp = head;
    do {
        std::cout << temp->data <<" "<< temp->spin << " ";
        temp = temp->next;
    } while (temp != head);
    std::cout << std::endl;
}
};



struct ContinuousTimeLattice{
    QfeLattice SpatialLattice; 
    UnionFind uf; 
    double invT; 
    std::vector<CircularLinkedList> TimeLadders; 
    std::vector<double> Tcouplings;  //For now, considering the simplified case where T couplings are uniform
    ContinuousTimeLattice(QfeLattice* spatialLattice, double inv_temp, double T_coup)
        : SpatialLattice(*spatialLattice), uf(), invT(inv_temp){
            Tcouplings.resize(SpatialLattice.n_sites);
            std::fill(Tcouplings.begin(), Tcouplings.end(), T_coup); 
            initializeTimeLadders(); 
        }

    private:
        void initializeTimeLadders(){
            int numLadders = SpatialLattice.n_sites;
            for (int i = 0; i< numLadders; i++){
                TimeLadders.emplace_back(invT); 
            } 
        }
};


bool checkWrap(Node* node){
    return node->data > node->next->data;
};
//Chatgpt
double intervalOverlap(double start1, double end1, double start2, double end2, double invT) {
    // Assume intervals are non-wrapping
    double overlapStart = std::max(start1, start2);
    double overlapEnd = std::min(end1, end2);
    return std::max(0.0, overlapEnd - overlapStart);
}
//Chatgpt
double returnOverlap(Node* node1, Node* node2, double invT){
    double nx = node1->data;
    double ny = node1->next->data;
    double mx = node2->data;
    double my = node2->next->data;

    double totalOverlap = 0.0;

    // Split intervals if they wrap
    std::vector<std::pair<double, double>> intervals1, intervals2;

    if (nx <= ny) {
        intervals1.push_back({nx, ny});
    } else {
        intervals1.push_back({nx, invT});
        intervals1.push_back({0.0, ny});
    }

    if (mx <= my) {
        intervals2.push_back({mx, my});
    } else {
        intervals2.push_back({mx, invT});
        intervals2.push_back({0.0, my});
    }

    // Compare all pairs
    for (auto& i1 : intervals1) {
        for (auto& i2 : intervals2) {
            totalOverlap += intervalOverlap(i1.first, i1.second, i2.first, i2.second, invT);
        }
    }

    return totalOverlap;
}


