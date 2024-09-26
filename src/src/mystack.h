#ifndef MYSTACK_H
#define MYSTACK_H

template<typename T>
class MyStack
{
private:
	T* storage;
	unsigned int capacity;
	unsigned int count;
	void ensureCapacity()
	{
		if (count == capacity)
		{
			capacity *= 2;
			T* temp = new T[capacity];
			for (int i = 0; i < count; i++)
			{
				temp[i] = storage[i];
			}
			delete storage;
			storage = temp;
		}
	}
public:
	MyStack()
	{
		capacity = 4;
		count = 0;
		storage = new T[capacity];
	}
	~MyStack()
	{
		delete[] storage;
	}
	void push(T val)
	{
		ensureCapacity();
		storage[count] = val;
		count++;
	}
	T pop()
	{
		count--;
		T ret = storage[count];
		return ret;
	}
	void clear()
	{
		count = 0;
	}
	bool isEmpty()
	{
		return count == 0;
	}
};

#endif // MYSTACK_H