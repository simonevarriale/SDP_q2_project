#ifndef BUCKETPQ_H
#define BUCKETPQ_H

#include <vector>
#include <map>
#include <algorithm>
#include <limits>

// Bucket Priority Queue
class BucketPriorityQueue {
private:
    std::map<int, std::vector<int>> buckets; // Map of priority to elements
    int minPriority; // Minimum priority value
    int maxPriority;

public:
    BucketPriorityQueue() {
        minPriority = std::numeric_limits<int>::max();
        maxPriority = std::numeric_limits<int>::min();
    }

    void print() {
        for (auto& entry : buckets) {
            auto& bucket = entry.second;
            std::cout << entry.first << ": ";
            for (auto it = bucket.begin(); it != bucket.end(); ++it) {
                std::cout << *it << " ";
            }
            std::cout << std::endl;
        }
    }

    // Retrieve the priority of a given element
    int getPriority(int element) const {
        for (const auto& entry : buckets) {
            if (std::find(entry.second.begin(), entry.second.end(), element) != entry.second.end()) {
                return entry.first;
            }
        }
        return -1; // Element not found in the priority queue
    }

    void deleteElement(int element) {
        auto priority = getPriority(element);
        if (priority == -1) {
            return;
        }
        auto& bucket = buckets[priority];
        for (int i = 0; i < bucket.size(); i++) {
            if (bucket[i] == element) {
                bucket.erase(bucket.begin() + i);
                if (bucket.empty()) {
                    buckets.erase(priority);
                    if (priority == maxPriority) {
                        maxPriority = std::numeric_limits<int>::min();
                        for (const auto& entry : buckets) {
                            maxPriority = std::max(maxPriority, entry.first);
                        }
                    }
                }
                break;
            }
        }
        // for (auto it = bucket.begin(); it != bucket.end(); it++) {
        //     if (*it == element) {

        //         it = bucket.erase(it);

        //         if (bucket.empty()) {
        //             buckets.erase(entry.first);
        //         }
        //         break;
        //     }
        // }
    }


    // Insert an element with a given priority
    void insert(int element, int priority) {
        if (buckets.find(priority) == buckets.end()) {
            buckets.insert({ priority, {} });
        }
        buckets[priority].push_back(element);
        minPriority = std::min(minPriority, priority);
        maxPriority = std::max(maxPriority, priority);
    }

    // Remove and return the element with the highest priority
    int extractMax() {
        int maxElement = buckets[maxPriority].back();
        buckets[maxPriority].pop_back();
        if (buckets[maxPriority].empty()) {
            buckets.erase(maxPriority);
            maxPriority = std::numeric_limits<int>::min();
            for (const auto& entry : buckets) {
                maxPriority = std::max(maxPriority, entry.first);
            }
        }
        return maxElement;
    }

    // Check if the bucket priority queue is empty
    bool isEmpty() const {
        return buckets.empty();
    }
};

#endif