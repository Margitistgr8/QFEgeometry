#pragma once

#include <cmath>
#include <set>
#include <stack>
#include <string>
#include <random>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <omp.h>       // OpenMP support
#include <mutex>       // If you use std::mutex for shared access
#include <utility>   
#include "CtLattice.h"


void PrintBinSegments(const std::unordered_map<int, std::vector<SpinSegment>>& bins,
                      int n_bins, double invT) {
    double bin_width = invT / n_bins;

    for (auto it = bins.begin(); it != bins.end(); ++it) {
        int bin_idx = it->first;
        double bin_start = bin_idx * bin_width;
        double bin_end = (bin_idx + 1) * bin_width;

        const std::vector<SpinSegment>& segs = it->second;

        std::cout << "Bin " << bin_idx << " ["
                  << bin_start << ", " << bin_end << ") contains "
                  << segs.size() << " segments:\n";

        for (size_t i = 0; i < segs.size(); ++i) {
            const SpinSegment& s = segs[i];
            std::cout << "  [" << s.start << ", " << s.end << "]  spin = " << s.spin << "\n";
        }
        std::cout << std::endl;
    }
}
void PrintBinSegmentsSim(const std::unordered_map<int, std::vector<SpinSegment>>& bins,
                      int n_bins, double invT) {
    double bin_width = invT / n_bins;

    for (auto it = bins.begin(); it != bins.end(); ++it) {
        int bin_idx = it->first;
        double bin_start = bin_idx * bin_width;
        double bin_end = (bin_idx + 1) * bin_width;

        const std::vector<SpinSegment>& segs = it->second;

        std::cout << "Bin " << bin_idx << " ["
                  << bin_start << ", " << bin_end << ") contains "
                  << segs.size() << " segments:\n";

        for (size_t i = 0; i < segs.size(); ++i) {
            const SpinSegment& s = segs[i];
            std::cout << "  [" << s.start << ", " << s.end << "]  spin = " << s.spin << "\n";
        }
        std::cout << std::endl;
    }
}
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

void DebugBuildRungSingleSite(int site_index) {
    CircularLinkedList& lad = lattice->TimeLadders[site_index];
    QfeSite& n1 = lattice->SpatialLattice.sites[site_index];
    double invT = lattice->invT;

    std::vector<SpinSegment> segs_i = extractSegments(lad, invT);
    std::cout << "Site " << site_index << " has " << segs_i.size() << " segments\n";

    for (int n2 = 0; n2 < n1.nn; ++n2) {
        int j = n1.neighbors[n2];
        CircularLinkedList& lad_j = lattice->TimeLadders[j];
        std::vector<SpinSegment> segs_j = extractSegments(lad_j, invT);

        double J = lattice->SpatialLattice.links[n1.links[n2]].wt;
        std::cout << "  Checking neighbor " << j << " (J = " << J << ")\n";

        for (const SpinSegment& s1 : segs_i) {
            for (const SpinSegment& s2 : segs_j) {
                if (s1.spin != s2.spin) continue;

                double overlap = intervalOverlap(s1.start, s1.end, s2.start, s2.end, invT);
                if (overlap > 0.0) {
                    double prob = 1.0 - exp(-2.0 * overlap * J);
                    std::cout << "    Overlap found: " << overlap 
                              << ", Spin = " << s1.spin 
                              << ", Prob = " << prob << "\n";
                }
            }
        }
    }
}
    

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
            }
            lad.mergeDuplicateSpins();
        }
    }




double Action(int N_THRESHOLD = 10){
    // Step 1: Prepare thread-local storage for deferred unions
    std::vector<std::pair<Node*, Node*>> all_union_pairs;
    
    //Use Binning to save time
    double avg_segment_length = 0.0;
    for (double gamma : lattice->Tcouplings) {
        avg_segment_length += 1.0 / gamma;
    }
    avg_segment_length /= lattice->Tcouplings.size();
    int n_bins = std::max(1, (int)(invT / avg_segment_length));

        // Global counters
    double globalAction = 0.0;
    // Step 2: Parallel outer loop over sites
     #pragma omp parallel
    {
        double localAction = 0.0; 
        #pragma omp for
        for (int i = 0; i < n_sites; i++) {
            CircularLinkedList& lad = lattice->TimeLadders[i];
            if (!lad.head) continue;


            QfeSite n1 = lattice->SpatialLattice.sites[i];
            std::vector<SpinSegment> segs_i = extractSegments(lad, invT);

            // Bin segments from site i
            std::unordered_map<int, std::vector<SpinSegment>> bins_i;
            bool bins_i_ready = false;

            for (int n2 = 0; n2 < n1.nn; n2++) {
                int j = n1.neighbors[n2];
                if (i >= j) continue;

                CircularLinkedList& nlad = lattice->TimeLadders[j];
                if (!nlad.head) continue;
                double wt = lattice->SpatialLattice.links[n1.links[n2]].wt;

                std::vector<SpinSegment> segs_j = extractSegments(nlad, invT);
                if (segs_i.size() <= N_THRESHOLD && segs_j.size() <= N_THRESHOLD) {
                    //printf("entering all to all step for spin site %d %d\n", i, j);
                    for (size_t s1 = 0; s1 < segs_i.size(); ++s1) {
                        for (size_t s2 = 0; s2 < segs_j.size(); ++s2) {
                            double overlap = intervalOverlap(segs_i[s1].start,segs_i[s1].end, segs_j[s2].start, segs_j[s2].end, invT);
                            if (overlap > 0.0) {localAction+= wt*overlap*segs_i[s1].spin*segs_j[s2].spin;}
                        }
                    }
                    continue;
                }
                //Build bins for site i if needed
                if (!bins_i_ready){
                for (size_t k = 0; k < segs_i.size(); ++k) 
                {
                std::vector<int> bins_spanned = getTimeBinsSpanned(segs_i[k].start, segs_i[k].end, invT, n_bins);
                for (int bin : bins_spanned) {bins_i[bin].push_back(segs_i[k]);}
                }
                bins_i_ready = true; 
                }

                //Build bins for site j if needed
                std::unordered_map<int, std::vector<SpinSegment>> bins_j;
                for (size_t k = 0; k < segs_j.size(); ++k) 
                {
                std::vector<int> bins_spanned = getTimeBinsSpanned(segs_j[k].start, segs_j[k].end, invT, n_bins);
                for (int bin : bins_spanned) {bins_j[bin].push_back(segs_j[k]);}
                }
                std::set<std::pair<Node*, Node*>> tried_pairs;                                
                for (int ind = 0; ind<n_bins; ind++) {
                    auto it1 = bins_i.find(ind);
                    auto it2 = bins_j.find(ind);
                    if (it1 == bins_i.end()||it2 == bins_j.end()) continue; 
                    double bin_width = invT / n_bins;
                    double bin_start = ind * bin_width;
                    double bin_end = (ind + 1) * bin_width;
                    const std::vector<SpinSegment>& list_i = it1->second;
                    const std::vector<SpinSegment>& list_j = it2->second;
                    for (size_t a = 0; a < list_i.size(); ++a) {
                        for (size_t b = 0; b < list_j.size(); ++b) {
                                const SpinSegment& s1 = list_i[a];
                                const SpinSegment& s2 = list_j[b];
                                 auto op = std::minmax(s1.node, s2.node);
                                if (tried_pairs.count(op)) continue; 
                                double overlap = intervalOverlap(s1.start, s1.end, s2.start, s2.end, invT);
                                if (overlap > 0.0) {
                                    tried_pairs.insert(op);
                                    localAction+=wt*overlap*s1.spin*s2.spin;
                                }
                        }
                    }

                }
            }
        }

        // Merge local unions into global list
        #pragma omp critical
        {
            globalAction += localAction;
        }
    }

    printf("Neighbor Action calculated %.12f \n", globalAction);
    //Adding transverse term
    for (int i=0; i<n_sites; i++)
     {
        CircularLinkedList& lad = lattice->TimeLadders[i]; 
        double gamma = lattice->Tcouplings[i];
        Node* temp = lad.head; 
        if (temp !=nullptr){
        do{
        double dt = lad.findInterval(temp);
        globalAction += (temp->spin)*dt*gamma; 
        temp = temp->next; 
        }while(temp!=lad.head);
        }
     }

    return globalAction;
}
/// @brief Create Cuts along the worldline
void CutLadder(std::vector<QfeRng>& thread_rngs){
    lattice->uf.reset();
    #pragma omp parallel
    {
        // Create thread-local RNG (clone of main one, or seeded separately)
        QfeRng& thread_rng = thread_rngs[omp_get_thread_num()];

        #pragma omp for
        for (int i=0; i<n_sites; i++)
        {
            CircularLinkedList& lad = lattice->TimeLadders[i]; 
            //printf("starting ladder %d with length %.12f\n", i, lad.findInterval(lad.head)); 
            if (lad.head == nullptr) continue; // Defensive: skip empty ladders

            Node* current = lad.head;
            double gamma = lattice->Tcouplings[i]; 
            std::vector<std::pair<Node*, std::vector<double>>> scheduledCuts;

            do{
                double intervalLength = lad.findInterval(current); 
                std::vector<double> cuts;
                std::vector<double> jumps; 
                double accumulated = 0.0;
                int count = 0; 
                do { 
                    double jump = rng.SampleExponential(gamma); 
                    //printf("starting %dst draw with length %.12f\n", count, jump); 
                    if (jump + accumulated > intervalLength){break;}
                    else{
                    accumulated+=jump; 
                    double cutTime = current->data + accumulated;
                    if (cutTime >= invT) cutTime -= invT;
                    cuts.push_back(cutTime);
                    jumps.push_back(jump);
                    count++; 
                    }
                }while(accumulated < intervalLength);
                if (!cuts.empty()) {
                    scheduledCuts.emplace_back(current, cuts);
                }

                current = current->next;
            }while(current!=lad.head);
            
            for (auto& pair : scheduledCuts) {
            Node* segment = pair.first;
            std::vector<double>& cuts = pair.second;
            Node* insertionPoint = segment;
            for (double cutTime : cuts) {
                lad.insertAfter(insertionPoint, cutTime);
                insertionPoint->next->spin = insertionPoint->spin;
                insertionPoint = insertionPoint->next;
                }
            }
        }
    }
// Add all nodes to union-find structure
    for (int i = 0; i < n_sites; i++) {
        CircularLinkedList& lad = lattice->TimeLadders[i];
        Node* temp = lad.head;
        if (!temp) continue;
        do {
            lattice->uf.addNode(temp);
            temp = temp->next;
        } while (temp != lad.head);
    }

}

std::vector<double> BuildRung(std::vector<QfeRng>& thread_rngs, int N_THRESHOLD=10){
    // Step 1: Prepare thread-local storage for deferred unions
    std::vector<std::pair<Node*, Node*>> all_union_pairs;
    double avg_segment_length = 0.0;
    for (double gamma : lattice->Tcouplings) {
        avg_segment_length += 1.0 / gamma;
    }
    avg_segment_length /= lattice->Tcouplings.size();
    int n_bins = std::max(1, (int)(invT / avg_segment_length));

        // Global counters
    int globalTotalAligned = 0;
    int globalConnectedAligned = 0;

    int globalcount = 0;
    double globalp = 0;
    double globaloverlap = 0;
    // Step 2: Parallel outer loop over sites
     #pragma omp parallel
    {
        QfeRng& thread_rng = thread_rngs[omp_get_thread_num()];        
        std::vector<std::pair<Node*, Node*>> local_unions;
        int localAligned = 0;
        int localconnectedAligned = 0;

        int     localcount = 0;
        double  localprob = 0.0;  
        double  localoverlap = 0.0; 


        #pragma omp for
        for (int i = 0; i < n_sites; i++) {
            CircularLinkedList& lad = lattice->TimeLadders[i];
            if (!lad.head) continue;


            QfeSite n1 = lattice->SpatialLattice.sites[i];
            std::vector<SpinSegment> segs_i = extractSegments(lad, invT);

            // Bin segments from site i
            std::unordered_map<int, std::vector<SpinSegment>> bins_i;
            bool bins_i_ready = false;

            for (int n2 = 0; n2 < n1.nn; n2++) {
                int j = n1.neighbors[n2];
                if (i >= j) continue;

                CircularLinkedList& nlad = lattice->TimeLadders[j];
                if (!nlad.head) continue;
                double wt = lattice->SpatialLattice.links[n1.links[n2]].wt;

                std::vector<SpinSegment> segs_j = extractSegments(nlad, invT);
                if (segs_i.size() <= N_THRESHOLD && segs_j.size() <= N_THRESHOLD) {
                    //printf("entering all to all step for spin site %d %d\n", i, j);
                    for (size_t s1 = 0; s1 < segs_i.size(); ++s1) {
                        for (size_t s2 = 0; s2 < segs_j.size(); ++s2) {
                            if (segs_i[s1].spin != segs_j[s2].spin) continue;
                                double overlap = intervalOverlap(segs_i[s1].start,segs_i[s1].end, segs_j[s2].start, segs_j[s2].end, invT);
                                if (overlap > 0.0) {
                                    double p = 1.0 - exp(-2.0 * overlap * wt);
                                    ++localcount;
                                    localprob += p;
                                    localoverlap+=overlap; 
                                    if (thread_rng.RandReal(0.0, 1.0) <= p) {
                                        local_unions.emplace_back(segs_i[s1].node, segs_j[s2].node);
                                        ++localconnectedAligned;
                                    }
                                }
                        }
                    }
                    continue; //double check if this is correct?
                }
                //Build bins for site i if needed
                if (!bins_i_ready){
                for (size_t k = 0; k < segs_i.size(); ++k) 
                {
                std::vector<int> bins_spanned = getTimeBinsSpanned(segs_i[k].start, segs_i[k].end, invT, n_bins);
                for (int bin : bins_spanned) {bins_i[bin].push_back(segs_i[k]);}
                }
                bins_i_ready = true; 
                }

                //Build bins for site j if needed
                std::unordered_map<int, std::vector<SpinSegment>> bins_j;
                for (size_t k = 0; k < segs_j.size(); ++k) 
                {
                std::vector<int> bins_spanned = getTimeBinsSpanned(segs_j[k].start, segs_j[k].end, invT, n_bins);
                for (int bin : bins_spanned) {bins_j[bin].push_back(segs_j[k]);}
                }
                std::set<std::pair<Node*, Node*>> tried_pairs;                                
                for (int ind = 0; ind<n_bins; ind++) {
                    //printf("entering binned step\n");
                    auto it1 = bins_i.find(ind);
                    auto it2 = bins_j.find(ind);
                    if (it1 == bins_i.end()||it2 == bins_j.end()) continue; 
                    const std::vector<SpinSegment>& list_i = it1->second;
                    const std::vector<SpinSegment>& list_j = it2->second;
                    for (size_t a = 0; a < list_i.size(); ++a) {
                        for (size_t b = 0; b < list_j.size(); ++b) {
                                const SpinSegment& s1 = list_i[a];
                                const SpinSegment& s2 = list_j[b];
                                 auto op = std::minmax(s1.node, s2.node);

                                if (s1.spin != s2.spin) continue;
                                if (tried_pairs.count(op)) continue; 
                                ++localAligned;

                                double overlap = intervalOverlap(s1.start, s1.end, s2.start, s2.end, invT);
                                if (overlap > 0.0) {
                                    tried_pairs.insert(op);
                                    double p = 1.0 - exp(-2.0 * overlap * wt);
                                    ++localcount;
                                    localprob += p;
                                    localoverlap+=overlap; 
                                    if (thread_rng.RandReal(0.0, 1.0) <= p) {
                                        local_unions.emplace_back(s1.node, s2.node);
                                        ++localconnectedAligned;
                                    }
                                }
                        }
                    }

                }
            }
        }
        #pragma omp atomic
        globalTotalAligned+=localAligned;

        #pragma omp atomic
        globalConnectedAligned+=localconnectedAligned;


         #pragma omp atomic
        globalcount+=localcount; 

        // Merge local unions into global list
        #pragma omp critical
        {
            globalp+=localprob;
            globaloverlap+=localoverlap;
            all_union_pairs.insert(all_union_pairs.end(),
                                   local_unions.begin(), local_unions.end());
        }
    }
    // Apply union operations serially
    for (size_t k = 0; k < all_union_pairs.size(); ++k) {
        lattice->uf.unionSets(all_union_pairs[k].first, all_union_pairs[k].second);
    }

    std::vector<double> res(3, -1.0);
    // if (globalTotalAligned != 0)
    // {res[0] = static_cast<double>(globalConnectedAligned) / globalTotalAligned;}
    if (globalcount != 0)
    {
    res[0] = static_cast<double>(globalConnectedAligned)/globalcount;
    res[1] = static_cast<double>(globalp) / globalcount;
    res[2] = static_cast<double>(globaloverlap) / globalcount;
    }
    return res;
}

void ExportWorldlinestoCSV(const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }
    file << "site_id,start_time,end_time,spin,cluster_id\n";

    for (int i = 0; i < n_sites; ++i) {
        CircularLinkedList& ladder = lattice->TimeLadders[i];
        Node* temp = ladder.head;

        if (!temp) continue;

        do {
            double start = temp->data;
            double end = temp->next->data;
            int spin = temp->spin;

            int cluster_id = -1;
            auto it = lattice->uf.nodeToId.find(temp);
            if (it != lattice->uf.nodeToId.end()) {
                cluster_id = lattice->uf.find(temp);
            }

            file << i << "," << start << "," << end << "," << spin << "," << cluster_id << "\n";

            temp = temp->next;
        } while (temp != ladder.head);
    }

    file.close();
    std::cout << "Worldline data exported to " << filename << std::endl;
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


int CountClusters() {
    std::unordered_set<int> unique_roots;

    for (int i = 0; i < n_sites; ++i) {
        CircularLinkedList& lad = lattice->TimeLadders[i];
        Node* temp = lad.head;
        if (!temp) continue;

        do {
            int root = lattice->uf.find(temp);
            unique_roots.insert(root);
            temp = temp->next;
        } while (temp != lad.head);
    }

    return static_cast<int>(unique_roots.size());
}

void UpdateSpins() {
    std::unordered_map<int, int> clusterSpin;

    // Step 1: assign a spin to each unique cluster root
    for (int i = 0; i < n_sites; i++) {
        CircularLinkedList& lad = lattice->TimeLadders[i];
        Node* temp = lad.head;
        if (!temp) continue;

        do {
            int root = lattice->uf.find(temp);
            if (clusterSpin.find(root) == clusterSpin.end()) {
                clusterSpin[root] = rng.RandBool() ? -1 : 1;
            }
            temp = temp->next;
        } while (temp != lad.head);
    }

    // Step 2: assign the spin to all nodes in each cluster
    for (int i = 0; i < n_sites; i++) {
        CircularLinkedList& lad = lattice->TimeLadders[i];
        Node* temp = lad.head;
        if (!temp) continue;

        do {
            int root = lattice->uf.find(temp);
            temp->spin = clusterSpin[root];
            temp = temp->next;
        } while (temp != lad.head);
    }

    lattice->uf.reset();

    // Step 3: clean up duplicate adjacent spins
    for (int i = 0; i < n_sites; i++) {
        CircularLinkedList& lad = lattice->TimeLadders[i];
        lad.mergeDuplicateSpins();
    }
}

int getSiteOfNode(Node* target) {
    for (int i = 0; i < n_sites; ++i) {
        CircularLinkedList& lad = lattice->TimeLadders[i];
        Node* temp = lad.head;
        if (!temp) continue;

        do {
            if (temp == target) return i;
            temp = temp->next;
        } while (temp != lad.head);
    }
    return -1;  // Not found
}

std::vector<double> BuildRungBug(std::vector<QfeRng>& thread_rngs, int N_THRESHOLD=0){
    // Step 1: Prepare thread-local storage for deferred unions
    std::vector<std::pair<Node*, Node*>> all_union_pairs;
    double avg_segment_length = 0.0;
    for (double gamma : lattice->Tcouplings) {
        avg_segment_length += 1.0 / gamma;
    }
    avg_segment_length /= lattice->Tcouplings.size();
    int n_bins = std::max(1, (int)(invT / avg_segment_length));

        // Global counters
    int globalTotalAligned = 0;
    int globalConnectedAligned = 0;

    int globalcount = 0;
    double globalp = 0;
    double globaloverlap = 0;
    QfeRng& thread_rng = thread_rngs[0];        
    std::vector<std::pair<Node*, Node*>> local_unions;
    int localAligned = 0;
    int localconnectedAligned = 0;

    int     localcount = 0;
    double  localprob = 0.0;  
    double  localoverlap = 0.0; 
    for (int i = 0; i < n_sites; i++) {
            CircularLinkedList& lad = lattice->TimeLadders[i];
            if (!lad.head) continue;


            QfeSite n1 = lattice->SpatialLattice.sites[i];
            std::vector<SpinSegment> segs_i = extractSegments(lad, invT);

            // Bin segments from site i
            std::unordered_map<int, std::vector<SpinSegment>> bins_i;
            bool bins_i_ready = false;

            for (int n2 = 0; n2 < n1.nn; n2++) {
                int j = n1.neighbors[n2];
                if (i >= j) continue;

                CircularLinkedList& nlad = lattice->TimeLadders[j];
                if (!nlad.head) continue;
                double wt = lattice->SpatialLattice.links[n1.links[n2]].wt;
            
                std::vector<SpinSegment> segs_j = extractSegments(nlad, invT);
                std::set<std::pair<Node*, Node*>> serial_pairs;
                std::set<std::pair<Node*, Node*>> binned_pairs;
                //  printf("ladder2 have size %d and %d\n", segs_i.size() , segs_j.size() );
                //printf("entering all to all step for spin site %d %d\n", i, j);
#if 0
                for (size_t s1 = 0; s1 < segs_i.size(); ++s1) {
                    for (size_t s2 = 0; s2 < segs_j.size(); ++s2) {
                        if (segs_i[s1].spin != segs_j[s2].spin) continue;
                            double overlap = intervalOverlap(segs_i[s1].start,segs_i[s1].end, segs_j[s2].start, segs_j[s2].end, invT);
                            if (overlap > 0.0) {
                                //printf("Spin Paring Spin Pairs of interval [%.12f, %.12f] and  [%.12f, %.12f]\n",segs_i[s1].start, segs_i[s1].end, segs_j[s2].start, segs_j[s2].end);
                                auto op = std::minmax(segs_i[s1].node, segs_j[s2].node);
                                serial_pairs.insert(op);
                                double p = 1.0 - exp(-2.0 * overlap * wt);
                                ++localcount;
                                localprob += p;
                                localoverlap+=overlap; 
                                if (thread_rng.RandReal(0.0, 1.0) <= p) {
                                    local_unions.emplace_back(segs_i[s1].node, segs_j[s2].node);
                                    ++localconnectedAligned;
                                }
                            }
                    }
                }
#endif
                //Build bins for site i if needed
                if (!bins_i_ready){
                for (size_t k = 0; k < segs_i.size(); ++k) 
                {
                    if (getSiteOfNode(segs_i[k].node)==-1){printf("Dangling Node is Found...\n");}
                std::vector<int> bins_spanned = getTimeBinsSpanned(segs_i[k].start, segs_i[k].end, invT, n_bins);
                for (int bin : bins_spanned) {bins_i[bin].push_back(segs_i[k]);}
                }
                bins_i_ready = true; 
                }

                //Build bins for site j if needed
                std::unordered_map<int, std::vector<SpinSegment>> bins_j;
                //printf("Bulding bins for companion ladder with size %d\n", segs_j.size());
                for (size_t k = 0; k < segs_j.size(); ++k) 
                {
                    if (getSiteOfNode(segs_j[k].node)==-1){printf("Dangling Node is Found...\n");}
                std::vector<int> bins_spanned = getTimeBinsSpanned(segs_j[k].start, segs_j[k].end, invT, n_bins);
                for (int bin : bins_spanned) {bins_j[bin].push_back(segs_j[k]);}
                }
                //PrintBinSegments(bins_j, n_bins, invT);
                std::set<std::pair<Node*, Node*>> tried_pairs;                                
                // printf("Considering Pair site %d %d\n", i,j);
                // printf("Printing Bin information for site %d:\n", i);
                // PrintBinSegments(bins_i, n_bins, invT);
                // printf("\n");
                // printf("Printing Bin information for site %d\n", j);
                // PrintBinSegments(bins_j, n_bins, invT);
                // printf("\n");
                for (int ind = 0; ind<n_bins; ind++) {
                    auto it1 = bins_i.find(ind);
                    auto it2 = bins_j.find(ind);
                    if (it1 == bins_i.end()||it2 == bins_j.end()) continue; 
                    double bin_width = invT / n_bins;
                    double bin_start = ind * bin_width;
                    double bin_end = (ind + 1) * bin_width;
                    const std::vector<SpinSegment>& list_i = it1->second;
                    const std::vector<SpinSegment>& list_j = it2->second;
                    //printf("Bin [%.12f, %.12f] scanned for rungs\n", bin_start, bin_end);

                    for (size_t a = 0; a < list_i.size(); ++a) {
                        for (size_t b = 0; b < list_j.size(); ++b) {


                                const SpinSegment& s1 = list_i[a];
                                const SpinSegment& s2 = list_j[b];
                        // std::cout<< " Considering interval: " ;
                        // printf(" [%.12f, %.12f] [%.12f, %.12f]  \n",s1.start,  s1.end, s2.start, s2.end);

                                 auto op = std::minmax(s1.node, s2.node);
                            //std::cout << "pointer 1 "<<s1.node << " pointer 2 "<<s2.node <<"\n";
                                if (s1.spin != s2.spin) continue;
                                if (tried_pairs.count(op)) 
                                {
                                    // printf("Previously counted\n");
                                    // std::cout << "Dumping all tried pairs (by address):\n";
                                    //     for (const auto& p : tried_pairs) {
                                    //         std::cout << "Pair: (" << p.first << ", " << p.second << ")\n";
                                    //     }
                                continue;} 
                                ++localAligned;

                                double overlap = intervalOverlap(s1.start, s1.end, s2.start, s2.end, invT);
                                if (overlap > 0.0) {
                                    //printf("Bin Paring Spin Pairs of interval [%.12f, %.12f] and  [%.12f, %.12f]\n", s1.start, s1.end, s2.start, s2.end);
                                    binned_pairs.insert(op);
                                    tried_pairs.insert(op);
                                    double p = 1.0 - exp(-2.0 * overlap * wt);
                                    ++localcount;
                                    localprob += p;
                                    localoverlap+=overlap; 
                                    if (thread_rng.RandReal(0.0, 1.0) <= p) {
                                        local_unions.emplace_back(s1.node, s2.node);
                                        ++localconnectedAligned;
                                    }
                                }
                        }
                    }

                }
            #if 0
            printf("testing all pairs in the serial bin---\n");
            int ii=0;
            for (auto& p: serial_pairs) {
                        ii++;
                        std::cout << ii << " site: (" << getSiteOfNode(p.first)
                                << ", " << getSiteOfNode(p.second)<<" ) ";
                        printf(" [%.12f, %.12f] [%.12f, %.12f]  \n",p.first->data,  p.first->next->data, p.second->data, p.second->next->data);

            }
            printf("testing all pairs in the binned bin---\n");
            ii=0;
            for (auto& p: binned_pairs) {
                        ii++;
                        std::cout << ii << " site: (" << getSiteOfNode(p.first)
                                << ", " << getSiteOfNode(p.second)<<" ) ";
                        printf(" [%.12f, %.12f] [%.12f, %.12f]  \n",p.first->data,  p.first->next->data, p.second->data, p.second->next->data);

            }

            if (serial_pairs != binned_pairs) 
            {
                std::cout << "Mismatch in pairings for sites " << i << " and " << j << "!\n";
                std::cout << "Serial only:\n";
                for (auto& p : serial_pairs) {
                    if (binned_pairs.count(p) == 0) {
                        std::cout << "  Missing in binned: (" << getSiteOfNode(p.first)
                                << ", " << getSiteOfNode(p.second) ;
                        //printf(" [%.12f, %.12f], [%.12f, %.12f]\n", p.first->data, p.first->next->data, p.second->data, p.second->next->data);
                        }    
                 }


                std::cout << "Binned only:\n";
                for (auto& p : binned_pairs) {
                    if (serial_pairs.count(p) == 0) {
                        std::cout << "  Extra in binned: (" << getSiteOfNode(p.first)
                                << ", " << getSiteOfNode(p.second);
                        //printf(" [%.12f, %.12f], [%.12f, %.12f]\n", p.first->data, p.first->next->data, p.second->data, p.second->next->data);
                        }
                }  
            }
            #endif



            }
        }

        globalTotalAligned+=localAligned;
        globalConnectedAligned+=localconnectedAligned;
        globalcount+=localcount; 
        {
            globalp+=localprob;
            globaloverlap+=localoverlap;
            all_union_pairs.insert(all_union_pairs.end(),
                                   local_unions.begin(), local_unions.end());
        }
    // Apply union operations serially
    for (size_t k = 0; k < all_union_pairs.size(); ++k) {
        lattice->uf.unionSets(all_union_pairs[k].first, all_union_pairs[k].second);
    }

    std::vector<double> res(3, -1.0);
    // if (globalTotalAligned != 0)
    // {res[0] = static_cast<double>(globalConnectedAligned) / globalTotalAligned;}
    if (globalcount != 0)
    {
    res[0] = static_cast<double>(globalConnectedAligned)/globalcount;
    res[1] = static_cast<double>(globalp) / globalcount;
    res[2] = static_cast<double>(globaloverlap) / globalcount;
    }
    return res;
    }


};