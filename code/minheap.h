#ifndef MINHEAP_H
#define MINHEAP_H

#define HEAP_FULL(n) (n == MAX_ELEMENTS-1)
#define HEAP_EMPTY(n) (!n)

extern int 	TotalNodes;
extern int	max_nodes;
extern Image	in_heap;

extern void ClearMinHeap();
extern void insert_min_heap(Node *item);
extern Node *delete_root();
extern void raise_to_root(int no);
extern void delete_item(int no);
extern Node *extract_min();
extern int In(Node *node);
extern int No_(Node *node);

#endif