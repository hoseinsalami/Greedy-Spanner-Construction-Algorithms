#ifndef SCAPEGOAT_H
#define SCAPEGOAT_H

template<class TKey, class TData>
struct ScapegoatNode {
	ScapegoatNode() {
		reset();
	}
	ScapegoatNode(TKey key, TData value) : key(key), value(value) {
		reset();
	}
	TKey key;
	TData value;
	ScapegoatNode<TKey, TData>* left;
	ScapegoatNode<TKey, TData>* right;
	ScapegoatNode<TKey, TData>* prev;
	ScapegoatNode<TKey, TData>* next;
	ScapegoatNode<TKey, TData>* parent;
	~ScapegoatNode() {
		if (left != nullptr)
			delete left;
		if (right != nullptr)
			delete right;
	}
	void reset() {
		left = nullptr;
		right = nullptr;
		prev = nullptr;
		next = nullptr;
		parent = nullptr;
	}
};

const double Alpha = 1.0 / 2.0;
const double Beta = 2;

template<class TKey, class TData>
class EmptyAugmentation {
public:
	static void process(ScapegoatNode<TKey, TData>* current) {
		/* Do nothing :) */
	}
};

template<class TKey, class TData, class TAugmentation = EmptyAugmentation<TKey, TData> >
class ScapegoatHeap {
private:
    unsigned int myCount, nrOfDeletedNodes;
    ScapegoatNode<TKey, TData>* root;
    ScapegoatNode<TKey, TData>* min;
    ScapegoatNode<TKey, TData>* max;
    unsigned int getSize(ScapegoatNode<TKey, TData>* node) {
        //return 1 + (Left != nullptr ? Left.GetSize() : 0) + (Right != nullptr ? Right.GetSize() : 0);
        ScapegoatNode<TKey, TData>* leftmost = node;
        while (leftmost->left != nullptr)
            leftmost = leftmost->left;
        ScapegoatNode<TKey, TData>* rightmost = node;
        while (rightmost->right != nullptr)
            rightmost = rightmost->right;
        unsigned int result = 1;
        while (leftmost != rightmost) {
            result++;
            leftmost = leftmost->next;
        }
        return result;
    }
    struct BalanceInfo {
        ScapegoatNode<TKey, TData>* top;
        ScapegoatNode<TKey, TData>* next;
    };
    void balance(ScapegoatNode<TKey, TData>* node, unsigned int n) {
        ScapegoatNode<TKey, TData>* leftmost = node;
        while (leftmost->left != nullptr)
            leftmost = leftmost->left;

        unsigned int i = 0;
        while ((((unsigned int)1) << i) - 1 < n)
            i++;
        //int leavesLeft = n - (int)Math.Pow(2, Math.Floor(Math.Log(n + 0.5, 2))) + 1;
        unsigned int leavesLeft = n - ((1 << (i - 1)) - 1);
        ScapegoatNode<TKey, TData>* currentParent = node->parent;
        BalanceInfo info;
        if (leavesLeft > 0) {
            leavesLeft--;
            info = balance(leftmost, 0, i - 1 /*(int)Math.Floor(Math.Log(n, 2))*/, leavesLeft);
        }
        else
            info = balance(leftmost, 1, i - 1 /*(int)Math.Floor(Math.Log(n, 2))*/, leavesLeft);

        if (currentParent == nullptr) {
            root = info.top;
            info.top->parent = nullptr;
        }
		else {
            info.top->parent = currentParent;
            if (currentParent->left == node)
                currentParent->left = info.top;
            else
                currentParent->right = info.top;
        }
		rebuildAugmentation(info.top);
		fixAugmentation(info.top->parent);
    }
    BalanceInfo balance(ScapegoatNode<TKey, TData>* node, unsigned int height, unsigned int maxHeight, unsigned int& leavesLeft) {
        nullChildren(node);
        ScapegoatNode<TKey, TData>* next;
        do {
            next = node->next;
            if (height == 1) {
                if (leavesLeft > 0) {
                    nullChildren(next);
                    addRightChild(node, next);
                    next = next->next;
                    leavesLeft--;
                }
                else
                    node->right = nullptr;
            }
            else if (height > 1) {
                BalanceInfo info;
                if (leavesLeft > 0) {
                    leavesLeft--;
                    info = balance(next, 0, height - 1, leavesLeft);
                }
                else
                    info = balance(next, 1, height - 1, leavesLeft);
                addRightChild(node, info.top);
                next = info.next;
            }
            if (height < maxHeight) {
                addParentLeft(node, next);
                node = next;
            }
            height++;
        }
        while (height <= maxHeight);

		BalanceInfo result;
		result.top = node;
		result.next = next;
        return result;
    }
    void nullChildren(ScapegoatNode<TKey, TData>* node) {
        node->left = nullptr;
        node->right = nullptr;
    }
    void addRightChild(ScapegoatNode<TKey, TData>* node, ScapegoatNode<TKey, TData>* child) {
        node->right = child;
        child->parent = node;
    }
    void addParentLeft(ScapegoatNode<TKey, TData>* node, ScapegoatNode<TKey, TData>* parent) {
        node->parent = parent;
        parent->left = node;
    }
	void finishInsert(ScapegoatNode<TKey, TData>* xin, ScapegoatNode<TKey, TData>* yin, ScapegoatNode<TKey, TData>* z, unsigned int depth)
    {
        ScapegoatNode<TKey, TData>* x = xin;
        ScapegoatNode<TKey, TData>* y = yin;
        if (y == nullptr) {
            root = z;
            min = z;
            max = z;
        }
        else {
            z->parent = y;
            if (z->key < y->key) {
                y->left = z;
                if (y->prev != nullptr) {
                    y->prev->next = z;
                    z->prev = y->prev;
                }
                else
                    min = z;
                y->prev = z;
                z->next = y;
            }
            else {
                y->right = z;
                if (y->next != nullptr) {
                    y->next->prev = z;
                    z->next = y->next;
                }
                else
					max = z;
                y->next = z;
                z->prev = y;
            }
        }
        if (depth > (log((double)myCount) / log(1.0 / Alpha))) {
            unsigned int height = 0;
            unsigned int leftSize = 0;
            unsigned int rightSize = 0;
            unsigned int size = 1;

            x = z;
            while (y != nullptr) {
                height++;
                if (y->left == x) {
                    leftSize = size;
                    rightSize = (y->right != nullptr) ? getSize(y->right) : 0;
                }
                else {
                    leftSize = (y->left != nullptr) ? getSize(y->left) : 0;
                    rightSize = size;
                }
                size = leftSize + rightSize + 1;

                if ((leftSize < Alpha * size) || (rightSize < Alpha * size)) {
                    balance(y, size);
                    break;
                }

                x = y;
                y = y->parent;
            }
        }
		else {
			fixAugmentation(z);
		}
    }
	void fixAugmentation(ScapegoatNode<TKey, TData>* node) {
		while (node != nullptr) {
			TAugmentation::process(node);
			node = node->parent;
		}
	}
	void rebuildAugmentation(ScapegoatNode<TKey, TData>* node) {
		if (node->left != nullptr)
			rebuildAugmentation(node->left);
		if (node->right != nullptr)
			rebuildAugmentation(node->right);
		TAugmentation::process(node);
	}
	inline void setParentRef(ScapegoatNode<TKey, TData>* z, ScapegoatNode<TKey, TData>* newValue) {
        if (z->parent->left == z)
            z->parent->left = newValue;
        else
            z->parent->right = newValue;
	}
	inline void checkRebalancingNeeded() {
		if (nrOfDeletedNodes * Beta > myCount) {
            if (myCount > 2)
                balance(root, myCount);
            nrOfDeletedNodes = 0;
        }
	}
public:
	ScapegoatHeap() {
		root = nullptr;
		min = nullptr;
		max = nullptr;
		myCount = 0;
		nrOfDeletedNodes = 0;
	}
	~ScapegoatHeap() {
		delete root;
	}
	template<class TResult>
    void augmentedSearch(TKey key, TResult& result) {
        ScapegoatNode<TKey, TData>* x = root;
        while (x != nullptr) {
			bool goingLeft = key <= x->key;
			result.update(x, goingLeft);
            if (goingLeft)
                x = x->left;
            else
                x = x->right;
        }
    }
    ScapegoatNode<TKey, TData>* search(TKey key) {
        ScapegoatNode<TKey, TData>* x = root;
        while (x != nullptr && key != x->key) {
            if (key < x->key)
                x = x->left;
            else
                x = x->right;
        }
        return x;
    }
    ScapegoatNode<TKey, TData>* searchClosest(TKey key) {
        ScapegoatNode<TKey, TData>* x = root;
        ScapegoatNode<TKey, TData>* p = root;
        while (x != nullptr && key != x->key) {
			p = x;
            if (key < x->key)
                x = x->left;
            else
                x = x->right;
        }
		if (x != nullptr) p = x;
        return p;
    }
    void remove(ScapegoatNode<TKey, TData>* z) {
        myCount--;
        nrOfDeletedNodes++;
		
		ScapegoatNode<TKey, TData>* x = nullptr;
		ScapegoatNode<TKey, TData>* y = nullptr;
        if (z->left == nullptr) {
            if (z->right == nullptr) {
                if (z->parent == nullptr) {
                    min = nullptr;
                    max = nullptr;
                    root = nullptr;
                    nrOfDeletedNodes = 0;
                    return;
                }
                else
                    setParentRef(z, nullptr);
            }
            else {
                z->right->parent = z->parent;
                if (z->parent != nullptr)
                    setParentRef(z, z->right);
                else
                    root = z->right;
            }
        }
        else {
            if (z->right == nullptr) {
                z->left->parent = z->parent;
                if (z->parent != nullptr)
                    setParentRef(z, z->left);
                else
                    root = z->left;
            }
            else {
				x = z->next;
				y = x->parent;
                if (x->right != nullptr) {
                    x->right->parent = x->parent;
                    setParentRef(x, x->right);
                }
                else
                    setParentRef(x, nullptr);

				z->next->prev = z->prev;
				z->prev->next = z->next;

				if (z == root)
					root = x;
				x->parent = z->parent;
				if (z->parent != nullptr)
                    setParentRef(z, x);
				x->left = z->left;
				if (z->left != nullptr)
					z->left->parent = x;
				x->right = z->right;
				if (z->right != nullptr)
					z->right->parent = x;

				checkRebalancingNeeded();
				fixAugmentation(y);
				fixAugmentation(x);
				return;
            }
        }

		if (z == min) {
			min = z->next;
			z->next->prev = nullptr;
		}
		else if (z == max) {
			max = z->prev;
			z->prev->next = nullptr;
		}
		else {
			z->next->prev = z->prev;
			z->prev->next = z->next;
		}

        checkRebalancingNeeded();
		fixAugmentation(z->parent);
    }
	unsigned int getCount() {
		return myCount;
	}
    ScapegoatNode<TKey, TData>* extractMin() {
        ScapegoatNode<TKey, TData>* result = min;
        remove(min);
        return result;
    }
    ScapegoatNode<TKey, TData>* getMin() {
        return min;
    }
    ScapegoatNode<TKey, TData>* extractMax() {
        ScapegoatNode<TKey, TData>* result = max;
        remove(max);
        return result;
    }
    ScapegoatNode<TKey, TData>* getMax() {
        return max;
    }
    ScapegoatNode<TKey, TData>* insert(TKey key, TData value) {
        myCount++;
        ScapegoatNode<TKey, TData>* z = new ScapegoatNode<TKey, TData>(key, value);
        ScapegoatNode<TKey, TData>* y = nullptr;
        ScapegoatNode<TKey, TData>* x = root;
        unsigned int depth = 0;
        while (x != nullptr) {
            y = x;
            if (key < x->key)
                x = x->left;
            else
                x = x->right;
            depth++;
        }
        finishInsert(x, y, z, depth);
		return z;
    }
    bool tryInsert(TKey key, TData& value) {
        myCount++;
        ScapegoatNode<TKey, TData>* z = new ScapegoatNode<TKey, TData>*(key, value);
        ScapegoatNode<TKey, TData>* y = nullptr;
        ScapegoatNode<TKey, TData>* x = root;
        unsigned int depth = 0;
        while (x != nullptr && key != x.Key) {
            y = x;
            if (key < x->key)
                x = x.Left;
            else
                x = x.Right;
            depth++;
        }
        if (x != nullptr && key == x.Key) {
            value = x->data;
            return false;
        }
        else {
            finishInsert(x, y, z, depth);
            return true;
        }
    }
};
#endif