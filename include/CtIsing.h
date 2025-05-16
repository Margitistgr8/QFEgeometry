#pragma once

#include <cmath>
#include <stack>
#include <string>
#include <random>
#include <unordered_map>
#include <vector>

#include "CtLattice.h"


class CtIsing {
    public: 
    ContinuousTimeLattice* lattice;
    double invT;     
    int n_sites;

    QfeRng rng; 


    CtIsing(ContinuousTimeLattice* lattice)
        : lattice(lattice), invT(lattice->invT), n_sites(lattice->SpatialLattice.n_sites) {} 


    void SeedRng(unsigned int seed) { rng = QfeRng(seed); }


    /// @brief Set all spins to +1
    void ColdStart() { 
        for (int i = 0; i<n_sites; i++)
        {   
            double start = rng.RandReal(0.0, invT);
            lattice->TimeLadders[i].insertAtEnd(start); 
            lattice->TimeLadders[i].head->spin = 1; 
        }
    };


    

    /// @brief Set all spins to random values
    void HotStart() { 
        for (int i = 0; i<n_sites; i++)
        {
            double counter = 0.0;
            CircularLinkedList& lad = lattice->TimeLadders[i]; 
            //Having a random starting point for each rung; 
            double start = rng.RandReal(0.0, invT); 
            lad.insertAtEnd(start); 
            lad.head->spin = rng.RandBool() ? -1 : 1;
            if (lad.head != nullptr)
            {
                //printf("Pass\n");
                Node* temp = lad.head;
                double gamma = lattice->Tcouplings[i]; 
                do
                {   double jump = rng.SampleExponential(gamma);
                    if (jump + counter > invT){
                        if (counter == 0.0)
                        {
                        temp->spin = rng.RandBool() ? -1 : 1;}
                        break;}
                    else{
                    counter+=jump;
                    double skip = start+counter< invT ? start+counter: start+counter-invT; 
                    lad.insertAtEnd(skip);  
                    temp = temp->next; 
                    temp->spin = rng.RandBool() ? -1 : 1;
                    }
                }while(counter < invT);
            // std::vector<double> ladder_intervals =  lad.returnAllIntervals();
            // double sum = std::accumulate(ladder_intervals.begin(), ladder_intervals.end(), 0.0);
            // double mean = (ladder_intervals.empty()) ? 0.0 : ladder_intervals.size()/sum;
            // printf("average segment length is: %.12f\n", mean);
            }
            lad.mergeDuplicateSpins();
        }
    }

    /// @brief Create Cuts along the worldline
    void CutLadder(){
        for (int i=0; i<n_sites; i++)
        {
            //printf("starting ladder %d\n", i); 
            CircularLinkedList& lad = lattice->TimeLadders[i]; 
            if (lad.head == nullptr) continue; // Defensive: skip empty ladders

            Node* current = lad.head;
            Node* head = lad.head; 
            double gamma = lattice->Tcouplings[i]; 
            do{
                double intervalLength = lad.findInterval(current); 
                double accumulated = 0.0; 
                Node* originalNext = current->next; // Save original next node (fixed target)

                Node* insertionPoint = current; // Start inserting from current node

                do { 
                    double jump = rng.SampleExponential(gamma); 
                    if (jump + accumulated > intervalLength){break;}
                    else{
                    accumulated+=jump; 
                    double startPosition = insertionPoint->data; 
                    double insertPosition = startPosition+jump< invT ? startPosition+jump: startPosition+jump-invT; 
                    lad.insertAfter(insertionPoint, insertPosition); 
                    insertionPoint->next->spin = insertionPoint->spin; 
                    insertionPoint = insertionPoint->next; 
                    }
                }while(accumulated < intervalLength);
                current = originalNext; 
            }while(current!=head);
        }
    }

void BuildRung(){

        //create an uf counter to avoid double counting. 
        AdjacencyMatrix tracker(n_sites); 
        
        //add all segments to the uf structure.
        for (int i=0; i<n_sites; i++)
        { 
            CircularLinkedList& lad = lattice->TimeLadders[i];
            Node* temp = lad.head; 
            do{
                lattice->uf.addNode(temp);
                temp = temp->next; 
            }while(temp!=lad.head); 
        }

        for (int i=0; i<n_sites; i++) //iterate over all spatial points
        {   
            //printf("Creating Rungs for site %d\n", i);
            CircularLinkedList& lad = lattice->TimeLadders[i]; 
            QfeSite n1 = lattice->SpatialLattice.sites[i]; 
            if (lad.head == nullptr) {continue;} // Skip this site if the list is empty
            for (int n2 = 0; n2<n1.nn; n2++) //iterate over neighbors
            {
                int j = n1.neighbors[n2];
                //printf("Looking at neighbor %d\n", j);
                int l = n1.links[n2]; 
                double wt = lattice->SpatialLattice.links[l].wt; 
                if (tracker.isConnected(i, j)){continue;}
                else{
                    
                    tracker.AddLink(i,j); 
                    CircularLinkedList& nlad = lattice->TimeLadders[j]; 
                    Node* temp = lad.head;
                    if (!temp) continue;
                    Node* neighbor = nlad.head;  
                    do{
                        //printf("checkpoint 1 at node %.12f \n", temp->data);
                        do
                        {   
                            double t_lap = returnOverlap(temp, neighbor, invT); 
                            if(t_lap> 0.0 && temp->spin==neighbor->spin){
                            //printf("checkpoint 2 with overlap %.12f and node %.12f\n", returnOverlap(temp, neighbor, invT), neighbor->data); 
                            double p = 1.0 - exp(-2.0*t_lap*wt); 
                            if (rng.RandReal(0.0, 1.0)< p)
                            {lattice->uf.unionSets(temp, neighbor);}}
                            neighbor = neighbor->next; 
                        } while (neighbor!=nlad.head);
                        temp=temp->next; 
                    }while(temp!=lad.head);  
                }
            }   
        }     
    }

double computeMeanSpin()
{
    double totalSpin = 0.0; 
     for (int i=0; i<n_sites; i++)
     {
        CircularLinkedList& lad = lattice->TimeLadders[i]; 
        Node* temp = lad.head; 
        if (temp !=nullptr){
        do{
        double dt = lad.findInterval(temp);
        totalSpin += (temp->spin)*dt; 
        temp = temp->next; 
        }while(temp!=lad.head);
        }
     }
    return totalSpin/(n_sites*invT);
}




void UpdateSpins(){
        std::vector<int> spins = rng.getRandomSpinArray(lattice->uf.nextId); 
        for (int i=0; i<n_sites; i++)
        { 
            CircularLinkedList& lad = lattice->TimeLadders[i];
            Node* temp = lad.head; 
            do{
                int label = lattice->uf.nodeToId[temp]; 
                temp->spin = spins[label]; 
                temp=temp->next; 
            }while(temp!=lad.head); 
        }

        lattice->uf.reset(); 

        //Deleting Potentially Duplicate Nodes
        for (int i=0; i<n_sites; i++)
        {
            CircularLinkedList& lad = lattice->TimeLadders[i]; 
            lad.mergeDuplicateSpins();
        }
    }
#if 0
    /// @brief Connecting worldline segments according probabilities
    void BuildRung(){

        //create an uf counter to avoid double counting. 
        AdjacencyMatrix tracker(n_sites); 
        
        //add all segments to the uf structure.
        for (int i=0; i<n_sites; i++)
        { 
            CircularLinkedList& lad = lattice->TimeLadders[i];
            Node* temp = lad.head; 
            do{
                lattice->uf.addNode(temp);
                temp = temp->next; 
            }while(temp!=lad.head); 
        }

        for (int i=0; i<n_sites; i++) //iterate over all spatial points
        {   
            //printf("Creating Rungs for site %d\n", i);
            CircularLinkedList& lad = lattice->TimeLadders[i]; 
            QfeSite n1 = lattice->SpatialLattice.sites[i]; 
            if (lad.head == nullptr) {continue;} // Skip this site if the list is empty
            for (int n2 = 0; n2<n1.nn; n2++) //iterate over neighbors
            {
                int j = n1.neighbors[n2];
                //printf("Looking at neighbor %d\n", j);
                int l = n1.links[n2]; 
                double wt = lattice->SpatialLattice.links[l].wt; 
                if (tracker.isConnected(i, j)){continue;}
                else{
                    tracker.AddLink(i,j); 
                    CircularLinkedList& nlad = lattice->TimeLadders[j]; 
                    Node* temp = lad.head;
                    Node* neighbor = nlad.head;
                    Node *nhead, *prevnode;  

                    Node* startNeighbor = neighbor;
                    do {
                        if (returnOverlap(temp, neighbor, invT) > 0.0) break;
                        neighbor = neighbor->next;
                    } while (neighbor != startNeighbor);

                    nhead = neighbor; 

                    do{
                        //printf("checkpoint 1 at node %.12f \n", temp->data);
                        Node* start = neighbor; 
                        while (returnOverlap(temp, neighbor, invT) > 0.0){
                        //printf("checkpoint 2 with overlap %.12f and node %.12f\n", returnOverlap(temp, neighbor, invT), neighbor->data); 
                        double p = 1.0 - exp(-2.0*t_lap*wt); 

                        if (rng.RandReal(0.0, 1.0)< p && temp->spin==neighbor->spin)
                        {lattice->uf.unionSets(temp, neighbor);}
                        prevnode = neighbor; 
                        neighbor = neighbor->next;
                        if (neighbor == start){break;}
                        }; 

                        //printf("checkpoint 3\n"); 
                        temp=temp->next; 
                        if (returnOverlap(temp, prevnode,invT)> 0.0){
                            //printf("checkpoint 4\n"); 
                            neighbor = prevnode;}
                    }while(temp!=lad.head);  

                    while(neighbor!=nhead)
                    {
                        if (rng.RandReal(0.0, 1.0)< .5 && temp->spin==neighbor->spin){lattice->uf.unionSets(temp, neighbor);}
                        neighbor = neighbor->next;
                    }; 
                }
            }   
        }     
    }

    /// @brief Update clusters of spins randomly
    void UpdateSpins(){
        std::vector<int> spins = rng.getRandomIntArray(lattice->uf.nextId); 
        for (int i=0; i<n_sites; i++)
        { 
            CircularLinkedList& lad = lattice->TimeLadders[i];
            Node* temp = lad.head; 
            do{
                int label = lattice->uf.nodeToId[temp]; 
                temp->spin = spins[label]; 
                temp=temp->next; 
            }while(temp!=lad.head); 
        }

        lattice->uf.reset(); 

        //Deleting Potentially Duplicate Nodes
        for (int i=0; i<n_sites; i++)
        {
            CircularLinkedList& lad = lattice->TimeLadders[i]; 
            Node* temp = lad.head;
            if (temp->next==lad.head){break;}
            do{     
              while(temp->next->spin == temp->spin)
              {
                lad.deleteNode(temp->next); 
              };
            temp = temp->next;
            }while(temp!=lad.head);
        }
    }
#endif
};

