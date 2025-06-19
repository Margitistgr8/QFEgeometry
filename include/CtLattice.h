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
            std::cout << "Node ID " << i << "with start, end: "<< node->data << " " << node->next->data <<" is in component with root ID " << root << std::endl;
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

    void   deleteNode(Node* node) {
    if (!node || !head) return;  // nothing to delete

    // Case 1: list has only one node
    if (head == node && head->next == head) {
        delete head;
        head = nullptr;
        member_num = 0;
        return;
    }

    // General case: find the node just before `node`
    Node* prev = head;
    bool found = false;

    do {
        if (prev->next == node) {
            found = true;
            break;
        }
        prev = prev->next;
    } while (prev != head);

    if (!found) {
        // Defensive: node not found in this list
        fprintf(stderr, "Warning: Attempted to delete node not in this list.\n");
        return;
    }

    // Rewire list: bypass the node to delete
    prev->next = node->next;

    // If deleting head, update it
    if (node == head) {
        head = node->next;
    }

    delete node;
    member_num--;
}
    double findInterval(Node* Node){//head 
        double val1  = Node->data;
        double val2  = Node->next->data;
        double dt = val2-val1; 

        if (dt<=0){return invT+dt;}
        else{return dt; } 
    }

    void mergeDuplicateSpins() {
    if (!head || head->next == head) return;

    Node* temp = head;
    do {
        // Stop if only one node remains
        if (temp == temp->next) break;

        // Delete duplicates while protecting circular structure
        while (temp->next != temp && temp->next->spin == temp->spin) {
            deleteNode(temp->next);  // relies on safe version of deleteNode
        }

        temp = temp->next;

    } while (temp != head);
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



//Chatgpt
std::vector<int> getTimeBinsSpanned(double start, double end, double invT, int n_bins) {
    std::vector<int> bins;
    double bin_width = invT / n_bins;

    auto timeToBin = [&](double t) {
        // clamp t to [0, invT)
        t = fmod(t + invT, invT);
        return static_cast<int>(t / bin_width);
    };

    if (start == end) {
        // Full circle: include all bins
        for (int b = 0; b < n_bins; ++b) {
            bins.push_back(b);
        }
    } else {
        int b_start = timeToBin(start);
        int b_end = timeToBin(end);

        if (end > start) {
            for (int b = b_start; b <= b_end; ++b) {
                bins.push_back(b % n_bins);
            }
        } else {
            // Wraparound: [start, invT) ∪ [0, end)
            for (int b = b_start; b < n_bins; ++b) {
                bins.push_back(b);
            }
            for (int b = 0; b <= b_end; ++b) {
                bins.push_back(b);
            }
        }
    }

    std::sort(bins.begin(), bins.end());
    bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
    return bins;
}

//Chatgpt
double intervalOverlap(double start1, double end1, double start2, double end2, double invT) {
    if (start1 == end1 && start2 == end2)
    return invT; //Both intervals completely wrap around worldline.
    
    //takes care of edge case when only one interval completely wraps around wordline
    if (start1 == end1) return (end2 >= start2) ? (end2 - start2) : (invT - start2 + end2);
    if (start2 == end2) return (end1 >= start1) ? (end1 - start1) : (invT - start1 + end1);


    auto segmentToRanges = [&](double start, double end) -> std::vector<std::pair<double, double>> {
        if (end >= start) {
            return {{start, end}};
        } else {
            return {{start, invT}, {0.0, end}};
        }
    };

    double overlap = 0.0;
    auto ranges1 = segmentToRanges(start1, end1);
    auto ranges2 = segmentToRanges(start2, end2);

    for (const auto& r1 : ranges1) {
        for (const auto& r2 : ranges2) {
            double s = std::max(r1.first, r2.first);
            double e = std::min(r1.second, r2.second);
            overlap += std::max(0.0, e - s);
        }
    }

    return overlap;
}
struct SpinSegment {
    double start;
    double end;
    int spin;
    Node* node;
};

void printSpinSegment(const SpinSegment& seg) {
    std::cout << "  Segment[start=" << seg.start
              << ", end=" << seg.end
              << ", spin=" << seg.spin
              << ", node_ptr=" << seg.node
              << "]\n";
}

std::vector<SpinSegment> extractSegments(CircularLinkedList& lad, double invT) {
    std::vector<SpinSegment> segments;
    if (!lad.head) return segments;

    Node* current = lad.head;
    do {
        double start = current->data;
        double end = current->next->data;
        //if (end < start) end += invT; // handle wraparound
        segments.push_back({start, end, current->spin, current});
        current = current->next;
    } while (current != lad.head);
    return segments;
}


