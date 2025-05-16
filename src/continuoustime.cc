#include <iostream>
#include "CtIsing.h"

using namespace std;

int main() {
    QfeLattice lattice; 
    lattice.InitRect(2,2, 1, 1); 
    ContinuousTimeLattice ctLattice(lattice, 10);
    std::fill(ctLattice.Tcouplings.begin(), ctLattice.Tcouplings.end(), 2);   
    CtIsing Ising(&ctLattice); 
    //Ising.ColdStart(); 
    Ising.HotStart(); 
    Ising.CutLadder(); 
            for (int i = 0; i<Ising.n_sites; i++)
        {
            printf("%d-----------\n", i);
            CircularLinkedList& lad = Ising.lattice->TimeLadders[i]; 
            Node* temp = lad.head; 
        if (temp != nullptr) {
        do {
            printf("%d %.12f\n", temp->spin, temp->data);
            temp = temp->next;
        } while (temp != lad.head && temp != nullptr);
        }
        printf("\n");
        }
    Ising.BuildRung();
    Ising.UpdateSpins();
    

    //  CircularLinkedList& lad = Ising.lattice->TimeLadders[0]; 
    //  CircularLinkedList& nlad = Ising.lattice->TimeLadders[2]; 
    //     Node* temp = lad.head;
    //     for (int i= 0; i<3; i++){temp= temp->next;}
    //     Node* neighbor = nlad.head;
    // cout<<"Checking Overlap "<< checkWrap(temp)<<endl;
    // cout<<"Checking Overlap "<< checkWrap(neighbor)<<endl; 
    // std::cout<<"overlap "<<returnOverlap(temp, neighbor, Ising.invT)<<endl;
    

    for (int i = 0; i<Ising.n_sites; i++)
        {
            printf("%d-----------\n", i);
            CircularLinkedList& lad = Ising.lattice->TimeLadders[i]; 
            Node* temp = lad.head; 
        if (temp != nullptr) {
        do {
            printf("%d %.12f\n", temp->spin, temp->data);
            temp = temp->next;
        } while (temp != lad.head && temp != nullptr);
        }
        printf("\n");
        }    
    // for (int i =0; i<16; i++){
    //     CircularLinkedList& lad =Ising.lattice->TimeLadders[i]; 
    //     Node* temp = lad.head; 
    //     if (temp != nullptr) {
    //     do {
    //         printf("%d ", temp->spin);
    //         temp = temp->next;
    //     } while (temp != lad.head && temp != nullptr);
    //     }
    //     printf("\n");
    // }
    // CircularLinkedList test = ctLattice.TimeLadders[0]; 
    // test.insertAtEnd(4); 
    // test.insertAtEnd(5); 
    // double val = test.findInterval(test.head);  
    // cout<<val<<"\n"<<endl; 
    // cout<<test.head->data<<endl; 
    return 0;
}

// int main(int argc, char *argv[])
// {

// return 0; 
// }

