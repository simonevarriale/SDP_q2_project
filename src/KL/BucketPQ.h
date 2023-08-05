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

    // Insert an element with a given priority
    void insert(int element, int priority) {
        buckets[priority].push_back(element);
        minPriority = std::min(minPriority, priority);
        maxPriority = std::max(maxPriority, priority);
    }

    // Remove and return the element with the highest priority
    int extractMax() {
        int maxElement = buckets[maxPriority].back();
        buckets.find(maxPriority)->second.pop_back();
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