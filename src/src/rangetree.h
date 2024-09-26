#ifndef RANGETREE_H
#define RANGETREE_H

#include <vector>

template<unsigned int dimP>
struct dpoint {
	int data;
	double coords[dimP];
	inline double get(unsigned int dimT) const {
		double result = coords[dimP - dimT];
		return result;
	}
	inline bool operator<=(const dpoint<dimP> &p) const {
		for (unsigned int i = 1; i <= dimP; i++) {
			if (get(i) > p.get(i)) {
				return false;
			}
		}
		return true;
	}
	inline double comp(const dpoint<dimP> &p, int dimT) const {
		if (get(dimT) != p.get(dimT))
			return get(dimT) < p.get(dimT);
		return data <= p.data;
	}
};

template<unsigned int dimP>
struct dimComp {
	unsigned int dimT;
	dimComp(unsigned int dimT) : dimT(dimT) {
	}
	inline bool operator()(const dpoint<dimP>& i, const dpoint<dimP>& j) {
		return i.comp(j, dimT); //i.get(dimT) < j.get(dimT);
	}
};

template<unsigned int dimP>
inline void pushbackIfAnswer(const dpoint<dimP>& point, const dpoint<dimP>& min, const dpoint<dimP>& max, std::vector<dpoint<dimP> >& result) {
	if (min <= point && point <= max) {
		result.push_back(point);
	}
}

template<unsigned int dimP>
class rangetree {
private:
	struct rangetreenode {
		unsigned int dimT;
		rangetreenode(std::vector<dpoint<dimP> >* points, unsigned int begin, unsigned int end) : dimT(1) {
			/*
			structure = nullptr;
			if (end - begin == 1) {
				point = points[begin];
				left = right = nullptr;
			}
			else {
				unsigned int middle = begin + (end + 1 - begin) / 2;
				point = points[middle - 1];
				left = new rangetreenode(points, begin, middle);
				right = new rangetreenode(points, middle, end);
			}
			*/
			this->points = points;
		}
		rangetreenode(std::vector<dpoint<dimP> >* points, std::vector<dpoint<dimP> >& pivots, unsigned int begin, unsigned int end) : dimT(2) {
			if (end - begin == 1) {
				point = (*points)[0];
				structure = nullptr;
				left = right = nullptr;
			}
			else {
				structure = new rangetree<dimP>(points, 1);
				unsigned int middle = begin + (end + 1 - begin) / 2;
				point = pivots[middle - 1];
				std::vector<dpoint<dimP> >* leftPoints = new std::vector<dpoint<dimP> >();
				std::vector<dpoint<dimP> >* rightPoints = new std::vector<dpoint<dimP> >();
				for (int i = 0; i < points->size(); i++) {
					if ((*points)[i].comp(point, 2))
						leftPoints->push_back((*points)[i]);
					else
						rightPoints->push_back((*points)[i]);
				}
				left = new rangetreenode(leftPoints, pivots, begin, middle);
				right = new rangetreenode(rightPoints, pivots, middle, end);
			}
		}
		rangetreenode(std::vector<dpoint<dimP> >* points, unsigned int begin, unsigned int end, unsigned int dimT) : dimT(dimT) {
			if (end - begin == 1) {
				point = (*points)[begin];
				structure = nullptr;
				left = right = nullptr;
			}
			else {
				std::vector<dpoint<dimP> >* newpoints = new std::vector<dpoint<dimP> >(points->begin() + begin, points->begin() + end);
				structure = new rangetree<dimP>(newpoints, dimT - 1);
				unsigned int middle = begin + (end + 1 - begin) / 2;
				point = (*points)[middle - 1];
				left = new rangetreenode(points, begin, middle, dimT);
				right = new rangetreenode(points, middle, end, dimT);
			}
		}
		~rangetreenode() {
			if (!isLeaf()) {
				delete structure;
				delete left;
				delete right;
			}
		}
		inline bool isLeaf() const {
			return left == nullptr;
		}
		int size() const {
			if (isLeaf()) return 1;
			return 1 + left->size() + right->size();
		}
		dpoint<dimP> point;
		rangetree<dimP>* structure;
		std::vector<dpoint<dimP> >* points;
		rangetreenode* left;
		rangetreenode* right;
		inline void reportAll(std::vector<dpoint<dimP> >& result) const {
			if (isLeaf())
				result.push_back(point);
			else {
				left->reportAll(result);
				right->reportAll(result);
			}
		}
		void query(const dpoint<dimP>& min, const dpoint<dimP>& max, std::vector<dpoint<dimP> >& result) {
			if (structure != nullptr)
				structure->query(min, max, result);
			else
				pushbackIfAnswer(point, min, max, result);
		}
		void queryAll(const dpoint<dimP>& min, const dpoint<dimP>& max, std::vector<dpoint<dimP> >& result) {
			if (isLeaf())
				pushbackIfAnswer(point, min, max, result);
			else {
				left->queryAll(min, max, result);
				right->queryAll(min, max, result);
			}
		}
		void query1D(const dpoint<dimP>& min, const dpoint<dimP>& max, std::vector<dpoint<dimP> >& result) {
			int minLow = 0;
			int minHigh = points->size()-1;
			while (minLow < minHigh) {
				int minMid = minLow + (minHigh - minLow)/2;
				if (min.get(1) <= (*points)[minMid].get(1))
					minHigh = minMid;
				else
					minLow = minMid + 1;
			}
			int maxLow = 0;
			int maxHigh = points->size()-1;
			while (maxLow < maxHigh) {
				int maxMid = maxLow + (maxHigh - maxLow)/2;
				if (max.get(1) < (*points)[maxMid].get(1))
					maxHigh = maxMid;
				else
					maxLow = maxMid + 1;
			}
			if (min.get(1) <= (*points)[minLow].get(1))
				result.push_back((*points)[minLow]);
			for (int i = minLow + 1; i < maxLow; i++)
				result.push_back((*points)[i]);
			if (max.get(1) >= (*points)[maxLow].get(1))
				result.push_back((*points)[maxLow]);
		}
	};
	rangetreenode* root;
	unsigned int dimT;
	void initialize(std::vector<dpoint<dimP> >* points, unsigned int dimT) {
		this->dimT = dimT;
		if (dimT == 1) {
			root = new rangetreenode(points, 0, points->size());
		}
		else if (dimT == 2) {
			std::vector<dpoint<dimP> > pivots(points->begin(), points->end());
			sort(pivots.begin(), pivots.end(), dimComp<dimP>(2));
			sort(points->begin(), points->end(), dimComp<dimP>(1));
			root = new rangetreenode(points, pivots, 0, points->size());
		}
		else {
			sort(points->begin(), points->end(), dimComp<dimP>(dimT));
			root = new rangetreenode(points, 0, points->size(), dimT);
		}
	}
	rangetreenode* findSplitNode(dpoint<dimP> min, dpoint<dimP> max) {
		rangetreenode* current = root;
		while (!current->isLeaf() && (max.get(dimT) <= current->point.get(dimT) || min.get(dimT) > current->point.get(dimT))) {
			if (max.get(dimT) <= current->point.get(dimT))
				current = current->left;
			else
				current = current->right;
		}
		return current;
	}
public:
	rangetree(std::vector<dpoint<dimP> >* points, unsigned int dimT) {
		initialize(points, dimT);
	}
	rangetree(std::vector<dpoint<dimP> >* points) {
		if (dimP == 1)
			sort(points->begin(), points->end(), dimComp<dimP>(1));
		initialize(points, dimP);
	}
	~rangetree() {
		delete root;
	}
	void query(const dpoint<dimP>& min, const dpoint<dimP>& max, std::vector<dpoint<dimP> >& result) {
	//	root->queryAll(min, max, result);
	//}
	//void query2(const dpoint<dimP>& min, const dpoint<dimP>& max, std::vector<dpoint<dimP> >& result) {
		if (dimT == 1)
			root->query1D(min, max, result);
		else {
			rangetreenode* splitNode = findSplitNode(min, max);
			if (splitNode->isLeaf())
				pushbackIfAnswer(splitNode->point, min, max, result);
			else {
				rangetreenode* current = splitNode->left;
				while (!current->isLeaf()) {
					if (min.get(dimT) <= current->point.get(dimT)) {
						current->right->query(min, max, result);
						current = current->left;
					}
					else {
						current = current->right;
					}
				}
				pushbackIfAnswer(current->point, min, max, result);

				current = splitNode->right;
				while (!current->isLeaf()) {
					if (max.get(dimT) >= current->point.get(dimT)) {
						current->left->query(min, max, result);
						current = current->right;
					}
					else {
						current = current->left;
					}
				}
				pushbackIfAnswer(current->point, min, max, result);
			}
		}
	}
};

#endif // RANGETREE_H
