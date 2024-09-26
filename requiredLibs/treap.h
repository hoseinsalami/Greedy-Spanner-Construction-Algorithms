#ifndef TREAP_H
#define TREAP_H

#include "random.h"

//static const unsigned int outofboundsindex = std::numeric_limits<unsigned int>::max();

template<class TKey, class TData>
struct Node {
	TData value;
	TKey key;
	Node<TKey, TData>* left;
	Node<TKey, TData>* right;
	Node<TKey, TData>* prev;
	Node<TKey, TData>* next;
	Node<TKey, TData>* parent;
	unsigned int priority;
	~Node() {
		if (left != NULL)
			delete left;
		if (right != NULL)
			delete right;
	}
};

template<class TKey, class TData>
class treap {
private:
	unsigned int count;
	Node<TKey, TData>* min;
	Node<TKey, TData>* max;
	Node<TKey, TData>* root;
	inline void balance(Node<TKey, TData>* index, Node<TKey, TData>* parent, unsigned int prio, bool isLeft) {
		while (parent != root) {
			unsigned int parentPrio = parent->priority;
			if (parentPrio < prio) {
				if (isLeft) {
					if (parent->right != NULL && parent->right->priority > prio) {
						index = parent->right;
						isLeft = false;
						prio = parent->right->priority;
					}
				}
				else {
					if (parent->left != NULL && parent->left->priority > prio) {
						index = parent->left;
						isLeft = true;
						prio = parent->left->priority;
					}
				}
				Node<TKey, TData>* parentParent = parent->parent;
				bool parentIsLeft = parentParent->left == parent;
				if (parentIsLeft)
					parentParent->left = index;
				else
					parentParent->right = index;
				
				parent->parent = index;
				index->parent = parentParent;
				Node<TKey, TData>* swappedTree;
				if (isLeft) {
					//rotate right
					swappedTree = index->right;
					parent->left = swappedTree;
					index->right = parent;
				}
				else {
					//rotate left
					swappedTree = index->left;
					parent->right = swappedTree;
					index->left = parent;
				}
				if (swappedTree != NULL)
					swappedTree->parent = parent;
				parent = parentParent;
				isLeft = parentIsLeft;
			}
			else
				return;
		}
		if (parent->priority < prio) {
			if (isLeft) {
				if (parent->right != NULL && parent->right->priority > prio) {
					index = parent->right;
					isLeft = false;
				}
			}
			else {
				if (parent->left != NULL && parent->left->priority > prio) {
					index = parent->left;
					isLeft = true;
				}
			}
			root = index;
			index->parent = NULL;
			parent->parent = index;
			Node<TKey, TData>* swappedTree;
			if (isLeft) {
				//rotate right
				swappedTree = index->right;
				parent->left = swappedTree;
				index->right = parent;
			}
			else {
				//rotate left
				swappedTree = index->left;
				parent->right = swappedTree;
				index->left = parent;
			}
			if (swappedTree != NULL)
				swappedTree->parent = parent;
		}
	}
public:
	inline treap() {
		count = 0;
		root = NULL;
	}
	inline Node<TKey, TData>* find(const TKey key) const {
		if (count == 0) return NULL;
		Node<TKey, TData>* index = root;
		while (index != NULL) {
			TKey currentKey = index->key;
			if (key == currentKey) return index;
			if (key < currentKey)
				index = index->left;
			else
				index = index->right;
		}
		return NULL;
	}
	inline Node<TKey, TData>* findClosest(const TKey key) const {
		Node<TKey, TData>* index = root;
		Node<TKey, TData>* prevIndex = NULL;
		while (index != NULL) {
			prevIndex = index;
			TKey currentKey = index->key;
			if (key == currentKey) return index;
			if (key < currentKey)
				index = index->left;
			else
				index = index->right;
		}
		return prevIndex;
	}
	inline void set(const TKey key, const TData value) {
		Node<TKey, TData>* index = root;
		Node<TKey, TData>* last = index;
		bool isLeft = false;
		while (index != NULL) {
			last = index;
			TKey currentKey = index->key;
			if (key == currentKey) {
				index->value = value;
				return;
			}
			isLeft = key < currentKey;
			if (isLeft)
				index = index->left;
			else
				index = index->right;
		}
		Node<TKey, TData>* oldRoot = root;
		if (count > 0) {
			Node<TKey, TData>* neu = new Node<TKey, TData>();
			neu->parent = last;
			neu->left = NULL;
			neu->right = NULL;
			if (isLeft) {
				neu->prev = last->prev;
				neu->next = last;
				last->prev = neu;
			}
			else {
				neu->prev = last;
				neu->next = last->prev;
				last->next = neu;
			}
			neu->key = key;
			neu->value = value;
			unsigned int prio = random_int();
			neu->priority = prio;
			if (isLeft)
				last->left = neu;
			else
				last->right = neu;
			count++;
			balance(neu, last, prio, isLeft);
		}
		else {
			count = 1;
			root = new Node<TKey, TData>();
			root->parent = NULL;
			root->left = NULL;
			root->right = NULL;
			root->prev = NULL;
			root->next = NULL;
			root->key = key;
			root->value = value;
			root->priority = random_int();
		}
	}
	inline void remove(const Node<TKey, TData>* index) {
		if (isValid(index->left)) {
			if (isValid(index->right)) {
				if (index->next != NULL) {

				}
				else {
				}
			}
			else {
				if (index->parent->left == index)
					index->parent->left = index->left;
				else
					index->parent->right = index->left;
				index->left->parent = index->parent;
				index->left = NULL;
			}
		}
		else {
			if (isValid(index->right)) {
				if (index->parent->left == index)
					index->parent->left = index->right;
				else
					index->parent->right = index->right;
				index->right->parent = index->parent;
				index->right = NULL;
			}
			else {
				if (index->parent->left == index)
					index->parent->left = NULL;
				else
					index->parent->right = NULL;
			}
		}
		if (index->prev != NULL)
			index->prev->next = index->next;
		if (index->next != NULL)
			index->next->prev = index->prev;
		delete index;
		count--;
	}
	inline void clear() {
		delete root;
		root = NULL;
		count = 0;
	}
	inline Node<TKey, TData>* getmin() {
		return min;
	}
	inline Node<TKey, TData>* getmax() {
		return max;
	}
};

#endif
