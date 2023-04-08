#ifndef _LINK_H_
#define _LINK_H_
#include "string"
  struct Link {
        int target, source;
        //std::string target_id,source_id;

        Link() {}

        Link(int s, int t) : source(s), target(t) {}

        friend bool operator<(const Link &n1, const Link &n2) {
            if (n1.source != n2.source) {
                return n1.source < n2.source;

            } else {
                return n1.target < n2.target;
            }

        }
    };
#endif