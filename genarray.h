#ifndef _GENARRAY_
#define _GENARRAY_

using namespace std;

template<typename T>
class genarray
{
	
	public:
		genarray();
		genarray(int newSize);
		genarray(const genarray& c);
		
		~genarray();
		
		int length();
		void resize(int newSize);
		
		genarray& operator= (const genarray& c);
		T& operator() (int k);
		T operator() (int k) const; 
		
	private:
		int myLength;
		T * myArray;
		
};


//function definitions:
template<typename T>
genarray<T>::genarray()
//default constructor - create a size 1 array.
{
	myLength = 1;
	myArray = new T[1];
}


template<typename T>
genarray<T>::genarray(int newSize)
//regular constructor
{
	myLength = newSize;
	myArray = new T[newSize];
}

template<typename T>
genarray<T>::genarray(const genarray& c)
//copy constructor
{
	myLength = c.myLength;
	myArray = new T[myLength];
	for (int i = 0; i < myLength; i++)
	{
		myArray[i] = c(i);
	}
}

template<typename T>
genarray<T>::~genarray()
//destructor
{
	delete [] myArray;
}

template<typename T>
int genarray<T>::length()
{
	return myLength;
}

template<typename T>
void genarray<T>::resize(int newSize)
{
	int numToCopy;
	if (myLength < newSize)
	{
		numToCopy = myLength;
	}
	else
	{
		numToCopy = newSize;
	}
	myLength = newSize;
	T * newArray;
	newArray = new T[newSize];
	for (int i = 0; i < numToCopy; i++)
	{
		newArray[i] = myArray[i];
	}
	myLength = newSize;
	delete [] myArray;
	myArray = newArray;
}

template<typename T>
genarray<T>& genarray<T>::operator= (const genarray& c)
//assignment operator
{
	myLength = c.myLength;
	T * newArray = new T[c.myLength];
	for (int i = 0; i < myLength; i++)
	{
		newArray[i] = c(i);
	}
	delete [] myArray;
	myArray = newArray;
}

template<typename T>
T& genarray<T>::operator() (int k)
{
	//if (k >= myLength)
	//{
		//cout << "genarray: out of bounds!" << endl;
		//cout << "value was: " << k << endl;
	//}
	//else
	//{
		return myArray[k];
	//}
}

template<typename T>
T genarray<T>::operator() (int k) const
{
	//if (k >= myLength)
	//{
		//cout << "genarray: out of bounds!" << endl;
		//cout << "value was: " << k << endl;
	//}
	//else
	//{
		return myArray[k];
	//}
}


#endif












	
