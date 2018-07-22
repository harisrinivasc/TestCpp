//
//  main.cpp
//  TestCPP
//
//  Created by Hari Chandrasekar on 6/21/18.
//  Copyright © 2018 Hari Chandrasekar. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <queue>
#include <stack>
#include <string>
#include <random>
//#include <cmath>
//#include "Sales_item.h"

//using std::string;
using namespace std;

//int reused = 42; // reused has global scope
typedef unsigned int uint;

class Book
{
    string Title;
    uint Price;
    uint NumSold;
public:
    Book() = default;
    Book(const string s, const uint price, const uint num) : Title(s), Price(price), NumSold(num) {};
    Book(const Book &copy)
    {
        Title = copy.Title;
        this->Price = copy.Price;
        this->NumSold = copy.NumSold;
    }
    void Change(const string s = "", const uint price = 99, const uint num = 99);
    friend istream& Read(istream&, Book&);
    //friend ostream& Print(ostream&, const Book&);
    string GetTitle(void) const;
    uint GetPrice(void) const;
    uint GetNumSold(void) const;
};

istream& Read(istream&, Book&);
ostream& Print(ostream&, const Book&);

istream& Read(istream& is, Book& Book1)
{
    cout << "Enter book Title, Price, NumSold " << endl;
    getline(is, Book1.Title);
    is >> Book1.Price >> Book1.NumSold;
    return is;
}

ostream& Print(ostream& os, const Book& Book1)
{
    os << "Title " << Book1.GetTitle() << ", Price " << Book1.GetPrice() << ", NumSold " << Book1.GetNumSold() << endl;
    return os;
}

void Book::Change(const string s, const uint price, const uint num)
{
    this->Title = s;
    this->NumSold = num;
    this->Price = price;
}

string Book::GetTitle(void) const
{
    return Title;
}

uint Book::GetPrice(void) const
{
    return Price;
}

uint Book::GetNumSold(void) const { return NumSold; }

int f(void);
int g(void);
int h(void);
int j(void);
bool BinSearch(const int, vector<int> *pvec);
int BinSearchRecursive(const int S, const vector<int> &A, const vector<int>::size_type lo, const vector<int>::size_type hi);
void WriteToQueue(queue<string> &WordQ, const string &Word);
long long int Factorial(const int N);
template <typename T> bool Popstack(stack<T> &S, T &Val);
void PrintSetlist(const list<unordered_set<int> > &setlist);
//bool operator< (const string &s1, const string &s2);
void PrintSortVec(const vector<int> &ivec);
void MyMerge(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
void MyMergeSort(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
void MyQuickSort(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
vector<int>::size_type MyQuickPartition(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
void SwapIvec(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi);
void ShuffleVec(vector<int> &ivec);
void HeapSwimUp(vector<int> &heaptree, vector<int>::size_type k);
void HeapSwimDown(vector<int> &heaptree, vector<int>::size_type k, vector<int>::size_type vecsize);
void InsertHeapTree(vector<int> &heaptree, int val);
int GetHeapTree(vector<int> &heaptree);

#define MAX_COL 5
#define MAX_ROW 3
//#define RC_TO_INT(R, C) (( ( R < MAX_ROW ? R : MAX_ROW - 1 ) * MAX_COL ) + ( C < MAX_COL ? C : MAX_COL - 1 ))
#define RC_TO_INT(R, C) ( ( ( R % MAX_ROW ) * MAX_COL ) + ( C % MAX_COL ) )

class classA
{
private:
    int temp;
public:
    classA()
    {
        temp = 0;
        cout << "Default constructor" << endl;
    }
    classA(const int temp_input) : temp(temp_input)
    {
        cout << "Constrcutor from int: " << temp_input << endl;
    }
    classA(const classA &class_input) : temp(class_input.temp)
    {
        cout << "Copy constructor: " << class_input.temp << endl;
    }
    void operator= (const classA& c1)
    {
        temp = c1.temp;
        cout << "Assignement: " << c1.temp << endl;
    }
    friend bool operator< (const classA& c1, const classA& c2);
    friend classA operator+ (const classA& c1, const classA& c2);
    friend bool operator== (const classA& c1, const classA& c2);
    friend bool operator!= (const classA& c1, const classA& c2);
    classA& operator+= (const classA& c1);
    virtual ~classA()
    {
        cout << "Destroying object: " << temp << endl;
    }
    int getTemp()
    {
        return temp;
    }
    void setTemp(int value)
    {
        temp = value;
    }
};

classA& classA::operator+= (const classA& c1)
{
    temp += c1.temp;
    return *this;
}
classA operator+ (const classA& c1, const classA& c2)
{
    classA sum;
    sum.temp = c1.temp + c2.temp;
    return sum;
}

bool operator< (const classA& c1, const classA& c2)
{
    return (c1.temp < c2.temp);
}

bool operator== (const classA& c1, const classA& c2)
{
    return (c1.temp == c2.temp);
}

bool operator!= (const classA& c1, const classA& c2)
{
    return !(c1 == c2);
}

/*
//MxN connected component matrix
#define GET_NUM(R,C)    (((R)*N)+(C))
const int M = 7, N = 6;
vector<pair<int,int>> MatrixNgbrs {{-1,0},{0,-1},{1,0},{0,1}};
vector<vector<int> > ComponentRoot(M,vector<int>(N,-1));
void SetRoot(const int R, const int C, const int RootVal);
const vector<vector<int>> ComponentInputs ={{0,1,1,0,0,0},
                                            {0,0,0,1,1,1},
                                            {1,1,1,0,0,0},
                                            {0,0,0,1,0,0},
                                            {0,1,0,0,0,1},
                                            {0,1,0,1,0,0},
                                            {0,1,1,1,0,0}};
*/

int NumWaysStair(const int Total, unordered_map<int,int> &ResultMap, const vector<int> &Steps);
int NumWaysNsum(const int Sum, const vector<int>::size_type VecIdx, unordered_map<string,int> &ResultMap, const vector<int> &InputVec);

class BstType;
class BstType
{
public:
	int key;
	string value;
	int size;
	shared_ptr<BstType> left;
	shared_ptr<BstType> right;
	BstType() : key(0), value(""), size(1), left(nullptr), right(nullptr) {}
	~BstType()
	{		
		cout << "Destruct: Key: " << key << " Value: " << value << " Size: " << size << endl;
	}
};

class GraphType
{
private:
	int NumVertex;
	int NumConnections;
	vector<unordered_set<int> > Connections;
public:
	GraphType(const int N)
	{
		NumVertex = N;
		NumConnections = 0;
		for(int i = 0; i < N; ++i)
			Connections.emplace_back(unordered_set<int> ());
	}
	unordered_set<int> &GetConnections(const int Vertex)
	{
		if (Vertex < NumVertex)
			return Connections[Vertex];
		else
			return (unordered_set<int>());
	}
	int GetDegree(const int Vertex)
	{
		if (Vertex < NumVertex)
			return Connections[Vertex].size();
		else
			return -9999;
	}
	void AddConnection(const int v1, const int v2)
	{
		if ((v1 < NumVertex) && (v2 < NumVertex))
		{
			//Add connection to both vertices
			auto ret1 = Connections[v1].insert(v2);
			auto ret2 = Connections[v2].insert(v1);
			//If either of the elements was not inserted (connection already present), do not
			//increment connection (edge) count
			if ((true == ret1.second) || (true == ret2.second))
				++NumConnections;
		}
	}
	void PrintConnections(void)
	{
		cout << "Print connections:" << endl;
		int cnt = 0;
		for (const auto &i : Connections)
		{
			cout << "Vertex: " << cnt++ << " connected to ";
			for (const auto &j : i)
				cout << j << " ";
			cout << endl;
		}
	}
};

class Dfs : public GraphType
{
private:
	vector<bool> Visited;
	vector<int> ParentVertex;
	const int Source;
	int count;
public:
	Dfs(const int N, const int SourceNode) : GraphType(N), Source(SourceNode)
	{
		for (int i = 0; i < N; ++i)
		{
			Visited.push_back(false);
			ParentVertex.push_back(i);
		}
	}
	void PerformDfs()
	{
		cout << "DFS trace from source: " << Source << endl;
		PerformDfs(Source);
		cout << endl;
	}
	void PerformDfs(const int v)
	{
		cout << v << " ";
		Visited[v] = true;
		++count;
		for(const auto &i : GetConnections(v))
			if (!Visited[i])
			{
				ParentVertex[i] = v;
				PerformDfs(i);
			}				
	}
};

class Bfs : public GraphType
{
private:
	vector<bool> Visited;
	vector<int> ParentVertex;
	const int Source;
	queue<int> VertexToVisit;
	int count;
public:
	Bfs(const int N, const int SourceNode) : GraphType(N), Source(SourceNode)
	{
		for (int i = 0; i < N; ++i)
		{
			Visited.push_back(false);
			ParentVertex.push_back(i);
		}
		VertexToVisit.push(Source);
	}
	void PerformBfs(void)
	{
		cout << "BFS trace from source: " << Source << endl;
		Visited[Source] = true;
		cout << Source << " ";
		while (!VertexToVisit.empty())
		{
			int v = VertexToVisit.front();
			VertexToVisit.pop();
			for (const auto &i : GetConnections(v))
				if (!Visited[i])
				{
					Visited[i] = true;
					ParentVertex[i] = v;
					++count;
					VertexToVisit.push(i);
					cout << i << " ";
				}
		}
		cout << endl;
	}
};

Dfs MyDfs(11, 0);
Bfs MyBfs(11, 0);

shared_ptr<BstType> Head = nullptr;
void InsertKey(const int key, const string value);
shared_ptr<BstType> InsertKey(shared_ptr<BstType> Node, const int key, const string value);
void DeleteKey(const int key);
shared_ptr<BstType> DeleteKey(shared_ptr<BstType> Node, const int key);
string GetValue(const int key);
void DeleteBst(void);
void DeleteBst(shared_ptr<BstType> Node);
void PrintBst(void);
void PrintBst(const shared_ptr<BstType> Node);
int GetSize(const shared_ptr<BstType> Node);
int GetMin(void);
shared_ptr<BstType> GetMin(const shared_ptr<BstType> Node);
int GetMax(void);
void DeleteMin(void);
shared_ptr<BstType> DeleteMin(shared_ptr<BstType> Node);
int GetRank(const int key);
int GetRank(shared_ptr<BstType> Node, const int key);

int main(int argc, const char * argv[])
{
    // insert code here...
    /*
    std::cout << "Hello, World!" << std::endl;
    std::cout << "Enter a number and a string: ";
    long int l_Test = 33;
    char s_Temp[20];
    std::cin >> l_Test >> s_Temp;
    std::cout << "Entered number and string: " << l_Test << " " << s_Temp << std::endl;
    int k = 0;
    for( int i = 0 ; i < 2 ; ++i, k++ )
        std::cout << i << " " << k << std::endl;
    std::cout << k << std::endl;
     */
    /*
    std::cout << "Enter numbers: " << std::endl;
    unsigned long int q_Input = 0, q_Sum = 0;
    if( std::cin >> q_Input )
    {
        while( std::cin >> q_Input )
        {
            std::cout << "Input number: " << q_Input << std::endl;
            q_Sum += q_Input;
        }
    }
    std::cout << "Sum :" << q_Sum << std::endl;
     */
    /*
    Sales_item book1, book2;
    cout << "Enter book1 and book 2: ";
    cin >> book1 >> book2;
    cout << "Book1 details " << book1 << endl;
    cout << "Book2 details " << book2 << endl;
    cout << "Book sum" << book1 + book2 << endl;
    while( cin >> book1 )
        cout << "Book loop " << book1 << endl;
    cerr << "Test error msg" << endl;
     */
    /*unsigned int q_A = 4, q_B = 10;
    int l_C = -40;
    cout << "A-B " << q_A - q_B << " B-A " << (int)q_B + l_C << endl;
    
    unsigned u = 10, u2 = 42;
    std::cout << u2 - u << std::endl;
    std::cout << u - u2 << std::endl;
    int i = 10, i2 = 42;
    std::cout << i2 - i << std::endl;
    std::cout << i - i2 << std::endl;
    std::cout << i - u << std::endl;
    
    char test_char[] = "A";
    cout << test_char << " test msg" " continue " << sizeof(test_char) << endl;
    */
    /*
    long double ld = 3.1415926536;
    int a(ld), b = (ld); // error: narrowing conversion required
    int c(ld), d = ld; // ok: but value will be truncated
    
    int units_sold1 = 0;
    int units_sold2 = {0};
    int units_sold3{0};
    int units_sold4(0);
    string test_string = "test string sdsfsffsdfsfsfsfsdfsdfsd";
    cout << test_string << " test msg" " continue " << sizeof(test_string) << endl;
    
    double x, y = 10.9;
    x = y = 10.9;
    auto xyz = pow(2,32) - 1;
    decltype(xyz) a2 = -4.7;
    cout << sizeof(xyz) << " " << sizeof(a2) << " " << xyz << " " << a2 << endl;
    
    double d1_img;
     */
    /*
    int unique = 0; // unique has block scope
    // output #1: uses global reused; prints 42 0
    std::cout << reused << " " << unique << std::endl;
    int reused = 0; // new, local object named reused hides global reused
    // output #2: uses local reused; prints 0 0
    std::cout << reused << " " << unique << std::endl;
    // output #3: explicitly requests the global reused; prints 42 0
    std::cout << ::reused << " " << unique << std::endl;
    */
    /*
    int i = 100, sum = 0;
    for (int i = 0; i != 10; ++i)
        sum += i;
    std::cout << i << " " << sum << std::endl;
     */
    /*
    int i = 1024, i2 = 2048; // i and i2 are both ints
    int &r = i, r2 = i2; // r is a reference bound to i; r2 is an int
    int i3 = 1024, &ri = i3; // i3 is an int; ri is a reference bound to i3
    int &r3 = r, &r4 = i2; // both r3 and r4 are references
    r3 = 10;
    cout << i << " " << sizeof(i) << " " << sizeof(r3) << endl;
    */
    /*
    int i = 1024, *p = &i, &r = i;
    int *&r1 = p;
    *r1 = 10;
    cout << i << endl;
    
    int ival = 1024;
    int *pi = &ival; // pi points to an int
    int **ppi = &pi; // ppi points to a pointer to an int
    */
    /*
    int i = 42;
    const int &r1 = i; // we can bind a const int& to a plain int object
    const int &r2 = 48; // ok: r1 is a reference to const
    const int &r3 = r1 * 2; // ok: r3 is a reference to const
    i++;
    cout << i << " " << r1 << " " << r2 << " " << r3 << endl;
    
    int a = 10;
    const int *ptr = &a;
    a++;
    cout << a << " " << *ptr << endl;
    */
    /*
    int i = -1, &r = 0, i2 = 2;
    int *const p2 = &i2;
    const int i5 = -1, &r5 = 0;
    const int *const p3 = &i2;
    const int *p1 = &i2;
    const int &const r2;
    const int i3 = i, &r4 = i;
    */
    /*
    int i = 0;
    int *const p1 = &i; // we can't change the value of p1; const is top-level
    const int ci = 42; // we cannot change ci; const is top-level
    const int *p2 = &ci; // we can change p2; const is low-level
    const int *const p3 = p2; // right-most const is top-level, left-most is not
    const int &r = ci; // const in reference types is always low-level
    
    i = ci; // ok: copying the value of ci; top-level const in ci is ignored
    p2 = p3; // ok: pointed-to type matches; top-level const in p3 is ignored
    
    int *p = p3; // error: p3 has a low-level const but p doesn't
    p2 = p3; // ok: p2 has the same low-level const qualification as p3
    p2 = &i; // ok: we can convert int* to const int*
    const int &r3 = ci; // error: can't bind an ordinary int& to a const int object
    const int &r2 = i; // ok: can bind const int& to plain int
    
    typedef double wages; // wages is a synonym for double
    wages base, *p; // base is a synonym for double, p for double*
    */
    /*
    int i = 0, &r = i;
    const int ci = i, &cr = ci;
    auto k = ci, &l = i; // k is int; l is int&
    auto &m = ci, *p = &ci; // m is a const int&;p is a pointer to const int
    // error: type deduced from i is int; type deduced from &ci is const int
    auto &n = i;
    auto *p2 = &ci;
    */
    /*
    // decltype of an expression can be a reference type
    int i = 42, *p = &i, &r = i;
    decltype(r + 0) b; // ok: addition yields an int; b is an (uninitialized) int
    decltype(*p) c; // error: c is int& and must be initialized
    */
    /*int i = 10;
    int &r = i;
    i++;
    r++;
    cout << i << " " << r << endl;
    
    constexpr const int *ptr = nullptr;
    (*ptr)++;
    */
    /*
    string str(10,'q'), str2;
    cout << str << " " << str.size() << " " << sizeof(str) << endl;
    //cin >> str >> str2;
    //cout << str << " " << str2 << " " << str.size() << " " << sizeof(str) << endl;
    
    //while( cin >> str )
        //cout << str << " " << str.size() << " " << sizeof(str) << endl;
    
    //while( getline(cin, str) )
    getline(cin, str);
    if( !str.empty() )
        cout << str << " " << str.size() << " " << sizeof(str) << endl;
    
    for( auto c : str )
        cout << c << " " << " " << sizeof(c) << endl;
    */
    /*
    string str;
    getline(cin, str);
    cout << str << " " << str.size() << " " << endl;
    string temp;
    for( auto c : str )
    {
        if ( !isspace(c) && isalnum(c) )
        {
            temp += c;
        }
        else
        {
            if( temp.size() )
            {
                cout << temp << " " << temp.size() << " " << endl;
            }
            temp = "";
        }
    }
    cout << temp << " " << temp.size() << " " << endl;
     */
    /*
    string s("some string");
    for (decltype(s.size()) index = 0; index != s.size() && !isspace(s[index]); ++index)
        s[index] = toupper(s[index]); // capitalize the current     character
    cout << s << endl;
    */
    /*
    string hex_num = "0123456789ABCDEF", output;
    string::size_type n;
    while( cin >> n )
    {
        if( n < hex_num.size() )
            output += hex_num[n];
    }
    cout << output << " " << output.size() << endl;
    */
    /*
    string c("test"), c2 = "test";
    string c3{"test"}, c4{'c','r', 'x'};
    cout << c << c2 << c3 << c4 << c4.size() << endl;
    */
    /*
    vector<string> v5{"hi"}; // list initialization: v5 has one element
    //vector<string> v6("hi"); // error: can't construct a vector from a string literal
    vector<string> v7{"1"}; // v7 has ten default-initialized elements
    vector<string> v8(10, "hi"); // v8 has ten elements with value "hi"
    cout << v5.size() << " " << v7.size() << " " << v8.size() << endl;
    cout << v5[0] << " " << v7[0] << " " << v8[0] << endl;
    return 0;
    
    int i = 10;
    auto &a = i;
     */
    /*
    string word;
    vector<string> svec;
    while (cin >> word)
    {
        svec.push_back(word);
    }
    //for( decltype(svec.size()) i = 0; i < svec.size(); ++i )
    for( auto i : svec )
        cout << i << " ";
    cout << endl << svec.size() << endl;
    */
    /*
    vector<int> v1 = {1,2,3,4}, v2 = {1,2,4,3};
    cout << (v1 == v2) << endl;
     */
    /*
    vector<int> input = {42, 65, 95, 100, 39, 67, 95, 76, 88, 76, 83, 92, 76, 93};
    vector<int> result(11,0);
    for( const auto i : input )
        if( i <= 100 )
            ++result[i/10];
    for( const auto i : result )
        cout << i << " ";
    cout << endl;
     */
    /*string s("test world");
    for( auto i = s.begin(); i != s.end() && !isspace(*i); ++i )
        *i = toupper(*i);
    cout << s << endl;
    
    int i = 10, *const ptr = &i;
    *ptr = 11;
    (*ptr)++;
    */
    /*
    vector<string> svec;
    string word;
    while( cin >> word )
        svec.push_back(word);
    for( const auto s : svec )
        //if( !s.empty() )
        if( isalnum( *s.cbegin() ) )
        {
            cout << s << " ";
            break;
        }
    cout << endl;
    */
    /*
    int i = 2 + 4 / 2;
    cout << i << endl;
     */
    /*
    unsigned cnt = 42; // not a constant expression
    constexpr unsigned sz = 42; // constant expression
    // constexpr see § 2.4.4 (p. 66)
    int arr[10]; // array of ten ints
    int *parr[sz]; // array of 42 pointers to int
    string bad[cnt]; // error: cnt is not a constant expression
    
    char a1[] = {'C', '+', '+'}; // list initialization, no null
    char a2[] = {'C', '+', '+', '\0'}; // list initialization, explicit null
    char a3[] = "C++"; // null terminator added     automatically
    //const char a4[6] = "Daniel"; // error: no space for the null!
    
    char ca[] = "C++\0";//{'C', '+', '+'}; // not null terminated
    cout << strlen(ca) << endl; // disaster: ca isn't null terminated
    */
    /*
    int i = 0, a[4][5];
    for( auto &row : a )
        for( auto &col : row )
            col = i++;
            
    for( const auto &row : a )
    {
        for( const auto col : row )
            cout << col << " ";
        cout << endl;
    }
    
    for( auto i = begin(a); i != end(a); ++i )
        for( auto j = begin(*i); j != end(*i); ++j )
            cout << *j << " ";
    cout << endl;
    
    for( auto i = begin(a); i != end(a); ++i )
        cout << *((int *)(*i) + 2) << " ";
    cout << endl;
    
    for( auto i = begin(a[3]); i != end(a[3]); ++i )
        cout << *i << " ";
    cout << endl;
    */
    /*
    int i;
    while((cin >> i) && (99 != i))
        cout << i << " ";
    cout << "exit loop" << endl;
     */
    /*
    typedef int int16;
    int16 a[10];
    using b = int[10];
    b c;
    cout << sizeof(b) << endl;
     */
    /*
    Book b1;
    Book b2("Harry Potter", 15, 10);
    b1.Change("Perry Mason", b2.GetPrice());
    Book b3;
    Book b4(b2);
    
    Read(cin, b3);
    
    Print(cout, b1);
    Print(cout, b2);
    cout << b3.GetTitle() << " " << b3.GetPrice() << " " << b3.GetNumSold() << endl;
    Print(cout, b4);
    
    cout << "End Test Class" << endl;
    */
    /*
    ifstream fin;
    string filename = "/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/TestInput.txt";
    fin.open(filename);
    vector<int> ivec;
    int temp;
    
    if(!fin)
        cerr << "Unable to open input file" << endl;
    else
    {
        cout << "file open success" << endl;
        while( fin >> temp )
            ivec.push_back(temp);
        fin.close();
        for(auto i : ivec)
            cout << i << endl;
    }
    
    ofstream fout;
    fout.open(filename + ".copy");
    if(!fout)
        cerr << "Unable to open output file" << endl;
    else
    {
        for(auto i = ivec.begin(); i != ivec.end(); ++i)
            fout << *i << endl;
    }
    fout.close();
     */
    /*
    //vector<int> ivec{2,4,5,6,8,10,14,15,19,20,30,40,50,60};
    vector<int> ivec{87,5,690,87,87,74,2,4,5,87,10,18,765,87,16,5,8,5,3,2,11};
    int search = 18;
    decltype(ivec.size()) e, b, m;
    bool elem_found = false;
    int cnt = 0;
    b = 0;
    e = ivec.size();
    m = (b+e)/2;
    while(b<e)
    {
        ++cnt;
        if(ivec[m] == search)
        {
            cout << "Element " << search << " found at index " << m << ". Cnt " << cnt << endl;
            elem_found = true;
            break;
        }
        else if(search < ivec[m])
            e = m;
        else
            b = m+1;
        m = (b+e)/2;
    }
    
    if( false == elem_found )
        cout << "Element " << search << " not found. Cnt " << cnt << endl;
    */
    /*
    vector<string> v = {"quasi", "simba", "frollo", "scar"};
    list<string> slist;
    // insert the last two elements of v at the beginning of slist
    cout << *slist.insert(slist.begin(), v.end() - 2, v.end()) << endl;
    slist.insert(slist.end(), {"these", "words", "will", "go", "at", "the", "end"});
    // run-time error: iterators denoting the range to copy from
    // must not refer to the same container as the one we are changing
    slist.insert(slist.begin(), slist.begin(), slist.end());
    */
    /*
    list<string> lst;
    string word;
    auto iter = lst.begin();
    while (cin >> word)
        iter = lst.insert(iter, word); // same as calling push_front
     */
    /*
    vector<int> ivec {1,2,3,4,5,6,7,8,9,10};
    auto iter = ivec.begin();
    while (iter != ivec.end())
    {
        if((*iter) % 2)
            iter = ivec.erase(iter);
        else
        {
            iter = ivec.insert(iter, *iter);
            iter += 2;
        }
    }
     */
    /*
    const char *cp = "Hello World!!!"; // null-terminated array
    char noNull[] = {'H', 'i'}; // not null terminated
    string s1(cp); // copy up to the null in cp; s1 == "Hello World!!!"
    string s2(noNull,2); // copy two characters from no_null; s2 == "Hi"
    string s3(noNull); // undefined: noNull not null terminated
    string s4(cp + 6, 5);// copy 5 characters starting at cp[6]; s4 == "World"
    string s5(s1, 6, 5); // copy 5 characters starting at s1[6]; s5 == "World"
    string s6(s1, 6); // copy from s1 [6] to end of s1; s6 == "World!!!"
    string s7(s1,6,20); // ok, copies only to end of s1; s7 == "World!!!"
    string s8(s1, 16); // throws an out_of_range exception
    */
    /*
    string s2 = "pi = 3.14";
    double d;
    // convert the first substring in s that starts with a digit, d = 3.14
    d = stod(s2.substr(s2.find_first_of("+-.0123456789")));
     */
    /*
    string line ("the quick red fox jumps over the slow red turtle back");
    vector<string> words;
    string::size_type b = 0;
    while( b < line.size() )
    {
        string::size_type space = line.find(" ",b);
        if(space != string::npos)
            words.push_back(line.substr(b, space-b));
        else
        {
            words.push_back(line.substr(b, line.size()-b));
            break;
        }
        b = space + 1;
    }*/
    
    /*
    vector<string> sort_words(words);
    sort(sort_words.begin(), sort_words.end());
    vector<string> sort_words_copy(sort_words);
    auto end_unique1 = unique(words.begin(), words.end());
    auto end_unique2 = unique(sort_words.begin(), sort_words.end());
    */
    
    /*
    //vector reverse
    vector<string> reverse_words(words);
    auto first = reverse_words.begin(), last = reverse_words.end();
    while((first!=last) && (first!=--last))
    {
        iter_swap(first,last);
        ++first;
    }
     */
    
    /*
    //string reverse
    string reverse_line(line);
    auto first = reverse_line.begin(), last = reverse_line.end();
    while((first!=last) && (first!=--last))
    {
        iter_swap(first,last);
        ++first;
    }
     */
    /*
    string line ("testline the quick/red fox*jumps over*the/slow red turtle back");
    //string line ("the quick red");
    string reverse_line;
    string::size_type end = line.size() - 1;
    while( (end > 0) && (end < line.size()) )
    {
        string::size_type space = line.find_last_of(" /",end);
        if(space != string::npos)
        {
            reverse_line += line.substr(space+1, end-space);
            //reverse_line += string(" ");
            reverse_line += line[space];
        }
        else
        {
            reverse_line += line.substr(0, end+1);
            break;
        }
        end = space - 1;
    }
    */
    /*
    vector<string> svec {"to", "another"};
    sort(svec.begin(), svec.end());
    */
    /*
    istream_iterator<int> in_iter(cin), eof; // read ints from cin
    vector<int> vec(in_iter, eof); // construct vec from an iterator range
    */
    
    /*
    // count the number of times each word occurs in the input
    map<string, int> word_count; // empty map from string to size_t
    string word;
    while (cin >> word)
        ++word_count[word]; // fetch and increment the counter for word
    for (const auto &w : word_count) // for each element in the map
        // print the results
        cout << w.first << " occurs " << w.second
        << ((w.second > 1) ? " times" : " time") << endl;
    */
    /*
    // count the number of times each word occurs in the input
    map<string, size_t> word_count; // empty map from string to size_t
    set<string> exclude = {"The", "But", "And", "Or", "An", "A",
        "the", "but", "and", "or", "an",
        "a"};
    string word;
    while (cin >> word)
        // count only words that are not in exclude
        if (exclude.find(word) == exclude.end())
            ++word_count[word]; // fetch and increment the counter for word
    
    for (const auto &w : word_count) // for each element in the map
        // print the results
        cout << w.first << " occurs " << w.second
        << ((w.second > 1) ? " times" : " time") << endl;
    */
    /*
    // define a vector with 20 elements, holding two copies of each number from 0 to 9
    vector<int> ivec;
    for (vector<int>::size_type i = 10; i != 0; --i) {
        ivec.push_back(static_cast<int>(i));
        ivec.push_back(static_cast<int>(i)); // duplicate copies of each number
    }
    // iset holds unique elements from ivec; miset holds all 20 elements
    set<int> iset(ivec.cbegin(), ivec.cend());
    multiset<int> miset(ivec.cbegin(), ivec.cend());
    cout << ivec.size() << endl; // prints 20
    cout << iset.size() << endl; // prints 10
    cout << miset.size() << endl; // prints 20
    */
    /*
    // count the number of times each word occurs in the input
    map<string, int> word_count; // empty map from string to size_t
    string word;
    while (cin >> word)
    {
        auto ret = word_count.insert(make_pair(word,1));
        if( false == ret.second )
            ++ret.first->second;
    }
    for (const auto &w : word_count) // for each element in the map
        // print the results
        cout << w.first << " occurs " << w.second
        << ((w.second > 1) ? " times" : " time") << endl;
    */
    /*
    map<string, int> word_count; // empty map from string to size_t
    cout << word_count["Anna"]; // fetch the element indexed by Anna; prints     1
    ++word_count["Anna"]; // fetch the element and add 1 to it
    cout << word_count["Anna"]; // fetch the element and print it; prints 2
    auto ret = word_count.find("Anna");
    cout << ret->second << (word_count.find("Anna"))->second;
    */
    /*
    ifstream rules("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/rules.txt");
    ifstream input("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/input.txt");
    ofstream output("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/output.txt");

    //generate rules map
    map<string, string> rules_map;
#if 0
    string rules_line;
    while(getline(rules, rules_line))
    {
        istringstream stream_from_line(rules_line);
        string key, value;
        while((stream_from_line >> key) && (stream_from_line >> value))
            rules_map[key] = value;
    }
#endif
    
    string key, value;
    while((rules >> key) && (rules >> value))
        rules_map[key] = value;
    
    //run through input file
    string input_line;
    while(getline(input, input_line))
    {
        istringstream stream_from_line(input_line);
        string word, new_word;
        while(stream_from_line >> word)
        {
            new_word = word;
            auto ret = rules_map.find(word);
            if(ret != rules_map.end())
                new_word = ret->second;
            output << new_word << " ";
        }
        output << endl;
    }
    
    rules.close();
    input.close();
    output.close();
        */
    /*
    // count the number of times each word occurs in the input
    unordered_map<string, int> word_count; // empty map from string to size_t
    string word;
    while (cin >> word)
        ++word_count[word]; // fetch and increment the counter for word
    for (const auto &w : word_count) // for each element in the map
        // print the results
        cout << w.first << " occurs " << w.second
        << ((w.second > 1) ? " times" : " time") << endl;
    */
    /*
    auto p = make_shared<ifstream>("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/rules.txt");
    p->close();
    */
    /*
    int *pi = new int(1024); // object to which pi points has value 1024
    string *ps = new string(10, '9'); // *ps is "9999999999"
    // vector with ten elements with values from 0 to 9
    vector<int> *pv = new vector<int>{0,1,2,3,4,5,6,7,8,9};
    */
    /*
    string word;
    vector<string> svec;
    
    cout << "enter words. type 'end' to terminate" << endl;
    while (cin>>word)
    {
        if (word == "end")
            break;
        svec.push_back(word);
    }
    
    if( !svec.empty() )
        for( auto &c : svec )
            c[0] = toupper(c[0]);
    
    for( const auto &s : svec )
        cout << s << " ";
    cout << endl;
    */
    /*
    const int sum = 8;
    vector<int> input_vec {1,2,3,5,8,9,10};
    cout << "enter num. type '-999' to terminate" << endl;
    
    int num;
    while(cin >> num)
    {
        if(-999 == num)
            break;
        input_vec.push_back(num);
    }
     */
    /*
    for (const auto i : input_vec)
        if( true == BinSearch(sum-i, &input_vec))
            cout << "element" << i << endl;
    */
    /*
    int ia[4][4] = { {0,1,2,3}, {4,5,6,7}, {8,9,10,11}, {12,13,14,15} };
    for( auto i:ia[2])
        cout << i << " ";
    cout << endl;
    for(auto i:ia)
        cout << *i << " ";
    cout << endl;
    for(auto i:ia)
        cout << *((int *)i + 2) << " ";
    cout << endl;
    for(int i = 0; i < 4; ++i)
        cout << **(ia + i) << " ";
    cout << endl;
    
    vector<int> ivec(begin(ia[0]), end(ia[0]));
    for(auto i:ivec)
        cout << i << " ";
    cout << endl;
    */
    /*
    cout << g() << f()+ (g()*h()) + j() << j() << endl;
    int i = 0;
    cout << i << " " << i++ << endl;
    
    short short_value = 32767;
    short_value += 1;
    cout << "short_value: " << short_value << endl;
    */
    /*
    ifstream inputf("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/input_para.txt");
    ofstream outputf("/Users/hari/C-programs/Xcode/TestCPP/TestCPP/TestCPP/output_para_srch_results.txt");
    
    string input_line, input_word;
    int linenum = 0, searchcnt = 0, wordnum = 0;
    //const string search_word = "memory";
    set<string> search_word = {"memory", "objects"};
    queue<string> WordQ;
    
    while (getline(inputf,input_line))
    {
        ++linenum;
        wordnum = 0;
        istringstream inputs_stream(input_line);
        while(inputs_stream >> input_word)
        {
            ++wordnum;
            WriteToQueue(WordQ, input_word);
            //trim input word to remove non-alphanumberic characters
            auto non_alnum = input_word.find_first_of(",.:;");
            if( string::npos != non_alnum )
                input_word.resize(non_alnum);
            if( search_word.end() != search_word.find(input_word) )
            {
                ++searchcnt;
                outputf << "  line(" << linenum << ") " << "word(" << wordnum << ") ";
                while(!WordQ.empty())
                {
                    outputf << WordQ.front() << " ";
                    WordQ.pop();
                }
                outputf << endl;
            }
        }
    }
    outputf << *search_word.cbegin() << "/" << *search_word.crbegin() << " occurs " << searchcnt << " times." << endl;
    
    inputf.close();
    outputf.close();
    */
    
    /*cout << RC_TO_INT(2,4) << endl;
    string sint("a1s2d4");
    cout << stoi(sint)+1 << endl;
    */
    /*
    classA c1(4), c2(c1), c4;
    c1.setTemp(9);
    classA c3 = c1 + c2;
    c3 += c1;
    c4 = c3;
    bool a = c1==c2;
    bool b = c1!=c2;
    bool c = c1<c2;
    cout << a << b << c << endl;
    */
    /*
    unordered_map<int,shared_ptr<classA>> class_map;
    int temp;
    while( cin >> temp )
    {
        auto temp_sptr = make_shared<classA>(temp);
        class_map[temp] = temp_sptr;
    }
    
    //for(const auto &i:class_map)
        //cout << i.first << " " << i.second << endl;
    for(auto i = class_map.cbegin(); i != class_map.cend(); ++i)
        cout << (*i).first << " " << (*i).second << endl;
    */
    /*
    //GCD
    int p,q,r=0,loopcnt=0;
    cin >> p >> q;
    while(0 != q)
    {
        ++loopcnt;
        r = p % q;
        p = q;
        q = r;
    }
    cout << "gcd: " << p << " loopcnt: " << loopcnt << endl;
     */
    /*
    vector<int> A;
    const int S = 88;
    const int N = 100;
    for(auto i = 0; i < N; ++i)
        A.push_back(2*i);
    auto hi = A.size();
    decltype(hi) lo = 0;
    cout << BinSearchRecursive(S, A, lo, hi) << endl;
     */
    /*
    int N;
    cin >> N;
    cout << Factorial(N) << endl;
    */
    
    /*
    //Djikstra's algm to solve expression with braces
    string expr = "(((2+(3*4))/2)";
    stack<int> oprand;
    stack<char> oprtr;
    
    for(const auto &i : expr)
    {
        if(isdigit(i))
            oprand.push(stoi(&i));
        else if(('*'==i) || ('/'==i) || ('+'==i) || ('-'==i))
            oprtr.push(i);
        else if(')'==i)
        {
            int a, b;
            char c;
            if(!Popstack(oprand,a))
            {
                cout << "Operand stack empty" << endl;
                break;
            }
            if(!Popstack(oprand,b))
            {
                cout << "Operand stack empty" << endl;
                break;
            }
            if(!Popstack(oprtr,c))
            {
                cout << "Operator stack empty" << endl;
                break;
            }
            switch(c)
            {
                case '+':
                {
                    oprand.push(a+b);
                    break;
                }
                case '-':
                {
                    oprand.push(b-a);
                    break;
                }
                case '*':
                {
                    oprand.push(a*b);
                    break;
                }
                case '/':
                {
                    oprand.push(b/a);
                    break;
                }
            }
        }
        else if('('==i)
        {
            //ignore
        }
        else
        {
            cout << "Invalid character: " << i << endl;
            break;
        }
    }
    
    int sol;
    bool ret = Popstack(oprand,sol);
    cout << "Solution: " << sol << " " << ret << endl;
    */
    
    /*
    //2SUM
    vector<int> ivec = {2,6,4,9,6,3,5,2,2,2,14,-3};
    int sum = 11;
    unordered_multiset<int> iset;
    
    for(const auto &i : ivec)
    {
        iset.insert(i);
        if(iset.end()!=iset.find(sum-i))
        {
            auto myrange = iset.equal_range(sum-i);
            while(myrange.first != myrange.second)
            {
                cout << "Pair found: " << i << " " << sum-i << endl;
                ++myrange.first;
            }
        }
    }
    */
    
    /*
    //3SUM
    vector<int> ivec = {2,6,4,9,6,3,5,2,2,2,14,-3};
    int sum = 11;
    unordered_multiset<int> iset;
    
    for(auto i = ivec.cbegin(); i != ivec.cend(); ++i)
    {
        iset.insert(*i);
        for(auto j = i+1; j < ivec.cend(); ++j)
        {
            int find_elem = sum-*i-*j;
            if(iset.end()!=iset.find(find_elem))
            {
                auto myrange = iset.equal_range(find_elem);
                while(myrange.first != myrange.second)
                {
                    cout << "3SUM Pair found: " << *i << " " << *j << " " << find_elem << endl;
                    ++myrange.first;
                }
            }
        }
    }
    */
    
    /*
    //Union find
    list<unordered_set<int> > setlist;
    const int N = 20;
    cout << "Max nodes: 1 to " << N << endl;
    cout << "1 Node#1 Node#2 for connecting" << endl;
    cout << "2 Node#1 Node#2 for checking" << endl;
    cout << "3 Print" << endl;
    int op = 0, node1 = 0, node2 = 0;
    while((cin >> op) &&
          ((3 == op) ||
           ((cin >> node1) && (cin >> node2))))
    {
        if((node2>N)||(node2>N)||(op<1)||(op>3))
        {
            cerr << "Invalid op: " << op << " Nodes: " << node1 << ", " << node2 << endl;
            break;
        }
            
        if(1==op)
        {
            auto node1itr = setlist.end(), node2itr = setlist.end();
            for(auto i = setlist.begin(); i != setlist.end(); ++i)
            {
                const unordered_set<int> &tempset = *i;
                
                if(tempset.end() != tempset.find(node1))
                    node1itr = i;
                if(tempset.end() != tempset.find(node2))
                    node2itr = i;
                if((node1itr != setlist.end()) && (node2itr != setlist.end()))
                    break;
            }
            if((node1itr != setlist.end()) && (node2itr != setlist.end()))
            {
                unordered_set<int> &tempset1 = *node1itr, &tempset2 = *node2itr;
                
                // Merge node2 set onto node1 set. Delete node2 set
                if(node1itr != node2itr)
                {
                    tempset1.insert(tempset2.cbegin(), tempset2.cend());
                    setlist.erase(node2itr);
                }
                else
                    cout << node1 << " and " << node2 << " are already connected" << endl;
            }
            else if((node1itr == setlist.cend()) && (node2itr != setlist.cend()))
            {
                unordered_set<int> &tempset2 = *node2itr;
                tempset2.insert(node1); //Add node1 to node2 set
            }
            else if((node1itr != setlist.cend()) && (node2itr == setlist.cend()))
            {
                unordered_set<int> &tempset1 = *node1itr;
                tempset1.insert(node2); //Add node2 to node1 set
            }
            else
            {
                //Node1 and node2 not present in setlist. Create new set and add to list
                unordered_set<int> tempset = {node1,node2};
                setlist.push_back(tempset);
            }
        }
        else if(2==op)
        {
            auto i = setlist.begin();
            for(; i != setlist.end(); ++i)
            {
                const unordered_set<int> &tempset = *i;
                if(tempset.end() != tempset.find(node1))
                {
                    //Found set with node1. Check if node2 is present in this set, and return
                    if(tempset.end() != tempset.find(node2))
                        cout << node1 << " and " << node2 << " are connected" << endl;
                    else
                        cout << node1 << " and " << node2 << " are NOT connected" << endl;
                    break;
                }
                if(tempset.end() != tempset.find(node2))
                {
                    //Found set with node2. Check if node1 is present in this set, and return
                    if(tempset.end() != tempset.find(node1))
                        cout << node1 << " and " << node2 << " are connected" << endl;
                    else
                        cout << node1 << " and " << node2 << " are NOT connected" << endl;
                    break;
                }
            }
            if(setlist.end() == i)
            {
                //Reached end of 'i' loop. This means Node1 and node2 were not found in the setlist
                cout << node1 << " and " << node2 << " are NOT connected" << endl;
            }
        }
        else if(3==op)
        {
            //Print setlist
            PrintSetlist(setlist);
        }
        
        //Clear op for while loop
        op = 0;
    }
    */
    
    /*
    //Selection sort
    //Find the min. element in the array and place it at the beginning
    //Const complexity of O(N) even if input is already sorted
    //typedef int SortType;
    //vector<SortType> ivec = {3,16,9,7,1,7,7,87,45,76};
    //vector<SortType> ivec = {"cat", "and", "zoo", "yes", "comic", "yanni", "able", "axle", "history"};
    //SortType temp;
    vector<int> ivec = {3,16,9,7,1,7,7,87,45,76};
    int temp;
    for(auto i = ivec.begin(); i != ivec.end(); ++i)
    {
        //Find min. in leftover array
        auto min = i;
        for(auto j = i; j != ivec.end(); ++j)
            if(*j < *min)
                min = j;
        //Swap minimum with LHS of array
        temp = *min;
        *min = *i;
        *i = temp;
    }
     
    for(const auto &i : ivec)
    cout << i << " ";
    cout << endl;
    */
    
    /*
    //Insertion sort
    //Traverse through each element and place it in the correct postion in the array
    //before it. If vector is already sorted, complexity is linear. No need to run 2nd loop
    //vector<int> ivec = {3,16,9,7,1,7,7,87,45,76};
    vector<int> ivec = {1,3,5,6,7,10,16,18,20,27,29};
    int temp;
    for(auto i = ivec.begin()+1; i != ivec.end(); ++i)
    {
        for(auto j = i; j != ivec.begin() && (*j < *(j-1)); --j)
        {
            temp = *j;
            *j = *(j-1);
            *(j-1) = temp;
            PrintSortVec(ivec);
        }
    }
     */
    
    /*
    //Merge sort
    //Sort 2 halves separately and then combine them. NlogN complexity
    //vector<int> ivec = {1,3,5,6,7,10,87,2,4,5,7,9,18,20,34};
    vector<int> ivec = {3,16,9,7,1,7,7,87,45,76};
    auto hi = ivec.size()-1;
    decltype(hi) lo = 0;
    MyMergeSort(ivec, lo, hi);
    */
    
    /*
    //Quick sort
    //Split/partition the array based on a partitioning element, and then move all items
    //less than the element to the left, and greater than the element to its right.
    //Sort the 2 partitions individually
    vector<int> ivec = {8,16,9,7,1,7,7,87,45,76};
    auto hi = ivec.size()-1;
    decltype(hi) lo = 0;
    MyQuickSort(ivec, lo, hi);
    */
    
    /*
    //Shuffle vector
    vector<int> ivec;
    for(int i = 0; i < 20; ++i)
    {
        ivec.push_back(i);
    }
    for(int i = 0; i < 10; ++i)
    {
        ShuffleVec(ivec);
        PrintSortVec(ivec);
    }
    */
    /*
    vector<int> ivec = {4,2,7,9,10,14,6,66,45,87,23,12,11,1,11,17,3};
    priority_queue<int> pq (ivec.cbegin(), ivec.cend());
    while(!pq.empty())
    {
        cout << pq.top() << " ";
        pq.pop();
    }
    cout << endl;
     */
    
    /*
    //Binary heap tree (max priority tree)
    vector<int> ivec = {4,2,7,9,10,14,6,66,45,87,23,12,11,1,11,17,3};
    ShuffleVec(ivec);
    PrintSortVec(ivec);
    vector<int> heaptree;
    for(const auto &i : ivec)
    {
        InsertHeapTree(heaptree, i);
        PrintSortVec(heaptree);
    }
    while(!heaptree.empty())
    {
        cout << "GetMax: " << GetHeapTree(heaptree) << " " << endl;
        PrintSortVec(heaptree);
    }
    */
    
    /*
    //Heapsort
    vector<int> ivec = {4,2,7,9,10,14,6,66,45,87,23,12,11,1,11,17,3};
    ShuffleVec(ivec);
    PrintSortVec(ivec);
    //construct heap upto half of the size
    for(auto i = ivec.size()/2; i < ivec.size(); --i)
        HeapSwimDown(ivec,i,ivec.size());
    PrintSortVec(ivec);
    
    //Swap Max element with end of array. Reduce array size and construct heap again
    auto vecsize = ivec.size();
    while(vecsize)
    {
        SwapIvec(ivec,0,vecsize-1);
        --vecsize;
        HeapSwimDown(ivec,0,vecsize);
    }
    PrintSortVec(ivec);
    */
    
    /*
    //Search for permutations of a string within another string
    const string s = "";
    const string y = "abcabckdakckkjbkbcaaabcbbc";
    unordered_multiset<char> srchset, srchset_copy;
    for(const auto &i : s)
        srchset.insert(i);
    srchset_copy = srchset;
    int srchcnt = 0;
    for(string::size_type i = 0; i < y.size(); ++i)
    {
        if(srchset_copy.find(y[i]) != srchset_copy.end())
        {
            ++srchcnt;
            auto range = srchset_copy.equal_range(y[i]);
            srchset_copy.erase(range.first); //Erase only 1 instance of character
        }
        else
        {
            srchcnt = 0;
            if(srchset_copy.size() != srchset.size())
                //Make a copy only if its not already the same
                srchset_copy = srchset;
        }
        if(srchcnt == s.size() && !s.empty())
        {
            cout << "Found: " << i << " " << y.substr(i+1-s.size(), s.size()) << endl;
            srchcnt = 0;
            srchset_copy = srchset;
            i += 1;
            i -= s.size(); //Move the loop counter backwards by s.size()-1 to find nested strings
        }
    }
    */
    
    /*
    //MxN connected component matrix - solution #1
#define GET_NUM(R,C)    (((R)*N)+(C))
    const int M = 7, N = 6;
    const int In[M][N] =   {{0,1,1,0,0,0},
                            {0,0,0,1,1,1},
                            {1,1,1,0,0,0},
                            {0,0,0,1,0,0},
                            {0,1,0,0,0,1},
                            {0,1,0,1,0,0},
                            {0,1,1,1,0,0}};
    list<unordered_set<int>> setlist;
    vector<pair<int,int>> Ngbrs {{-1,0},{0,-1},{1,0},{0,1}};
    
    //Add to setlist
    for(int i = 0; i < M; ++i)
        for(int j = 0; j < N; ++j)
            if(In[i][j])
            {
                const int CurrNode = GET_NUM(i,j);
                vector<decltype(setlist.begin())> NgbrItr;
                for(auto &k : Ngbrs)
                {
                    int NgbrRow = i+k.first, NgbrCol = j+k.second;
                    if((NgbrRow<0) || (NgbrRow>=M)||
                       (NgbrCol<0) || (NgbrCol>=N))
                        continue;
                    int Ngbr = GET_NUM(NgbrRow, NgbrCol);
                    for(auto l = setlist.begin(); l != setlist.end(); ++l)
                    {
                        unordered_set<int> &TempSet = *l;
                        if(TempSet.find(Ngbr) != TempSet.end())
                            NgbrItr.push_back(l);
                    }
                }
                
                if(NgbrItr.empty())
                {
                    //Input not found. Create a new set (component) and add to the list
                    unordered_set<int> TempSet = {CurrNode};
                    setlist.push_back(TempSet);
                    continue;
                }
                if(NgbrItr.size() > 1)
                {
                    //More than 1 ngbr is already present in the setlist. Merge all of them into 1 set and
                    //delete the set that was merged
                    for(auto l = 1; l < NgbrItr.size(); ++l)
                    {
                        NgbrItr[0]->insert(NgbrItr[l]->begin(), NgbrItr[l]->end());
                        setlist.erase(NgbrItr[l]);
                    }
                }
                
                //Add the new element to the setlist
                NgbrItr[0]->insert(CurrNode);
            }
    */
    
    /*
    //MxN connected component matrix - solution #2
    int ComponentCnt = 0;
    for(int i = 0; i < M; ++i)
        for(int j = 0; j < N; ++j)
            if((ComponentInputs[i][j]) && (-1 == ComponentRoot[i][j]))
            {
                ++ComponentCnt;
                SetRoot(i,j, ComponentCnt);
            }
    
    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < N; ++j)
            cout << ComponentRoot[i][j] << " ";
        cout << endl;
    }
    */
    
    /*
    //Reverse vowels of a string
    string s = "xcvtr";
    vector<pair<char,int>> VowelPos;
    const unordered_set<char> Vowels = {'a','e','i','o','u'};
    for(string::size_type i = 0; i < s.size(); ++i)
        if(Vowels.find(s[i]) != Vowels.end())
            VowelPos.push_back(make_pair(s[i], i));
    for(decltype(VowelPos.size()) i = 0; i < VowelPos.size(); ++i)
        s[VowelPos[i].second] = VowelPos[VowelPos.size()-1-i].first;
    */
    
    /*
    //Atom counter
    const string formula = "(H2O2He3Mg4O4)2";
    map<string,int> AtomMap;
    string AtomName;
    int AtomCnt = 0;
    bool IsParanthesis = false;
    for(const auto &i : formula)
    {
        if(isalpha(i))
        {
            if(isupper(i))
            {
                //Start New atom. Before that, store previous atom in map
                if(!AtomName.empty())
                {
                    auto MapItr = AtomMap.find(AtomName);
                    if(AtomMap.find(AtomName) == AtomMap.end())
                        //Atom not present. Create new element
                        AtomMap[AtomName] = AtomCnt;
                    else
                        //Atom already exists. Add to previous cnt
                        MapItr->second += AtomCnt;
                }
                AtomName = i;
                AtomCnt = 1;
            }
            else
                //Atom name has multiple characters
                AtomName += i;
        }
        if(isdigit(i))
        {
            int Num = static_cast<int>(i) - static_cast<int>('0');
            if(IsParanthesis)
            {
                //Multiply all atom counts in map
                for(auto &j : AtomMap)
                    j.second *= Num;
                IsParanthesis = false;
            }
            else
                AtomCnt = Num;
        }
        if(i == ')')
        {
            IsParanthesis = true;
            //Store previous atom in map
            if(!AtomName.empty())
            {
                auto MapItr = AtomMap.find(AtomName);
                if(AtomMap.find(AtomName) == AtomMap.end())
                    //Atom not present. Create new element
                    AtomMap[AtomName] = AtomCnt;
                else
                    //Atom already exists. Add to previous cnt
                    MapItr->second += AtomCnt;
            }
        }
    }
    */
    
    /*
    //Staircase problem (coin problem)
    const int N = 7; //need to reach N steps using any combination of below steps
    vector<int> Steps = {5,3,1};
    unordered_map<int,int> ResultMap; //Temporary space to store results to avoid multiple recursion
    cout << NumWaysStair(N,ResultMap,Steps) << endl;
    for(const auto &i : ResultMap)
        cout << i.first << " " << i.second << endl;
    */
    
    /*
    //Find the Longest Increasing Subsequence
    //const vector<int> InputVec = {10,8,3,9,2,11,4,7,12,13,14};
    const vector<int> InputVec = {6,4,5,3,2,1,7,8,9,10,11};
    list<set<int>> SetList;
    for(const auto &i : InputVec)
    {
        bool InputAdded = false;
        set<int>::size_type PrevSetSize = 0;
        for(auto j = SetList.begin(); j != SetList.end(); ++j)
        {
            set<int> &TempSet = *j;
            if(i >= *TempSet.rbegin())
            {
                if(TempSet.size() == PrevSetSize)
                    //Delete this set as the same element is added to 2 sets of same size
                    SetList.erase(j);
                else
                {
                    TempSet.insert(i);
                    InputAdded = true;
                    PrevSetSize = TempSet.size() - 1;
                }
            }
        }
        if(!InputAdded)
            //Create a new set
            SetList.emplace_back(set<int>{i});
    }
    set<int>::size_type MaxSize = 0;
    for(const auto &i : SetList)
    {
        const set<int> &TempSet = i;
        if(MaxSize < TempSet.size())
            MaxSize = TempSet.size();
    }
    for(const auto &i : SetList)
    {
        const set<int> &TempSet = i;
        if(MaxSize == TempSet.size())
        {
            for(const auto &j : TempSet)
                cout << j << " ";
            cout << endl;
        }
    }
    */
    
    /*
    //Num ways to get sum (N-sum)
    const int Sum = 16; //need to reach N steps using any combination of below steps
    vector<int> InputVec = {6,4,10,2};
    unordered_map<string,int> ResultMap; //Temporary space to store results to avoid multiple recursion
    cout << NumWaysNsum(Sum,InputVec.size()-1,ResultMap,InputVec) << endl;
    for(const auto &i : ResultMap)
        cout << i.first << " " << i.second << endl;
    */
    
    	/*
	//Binary search tree BST
	cout << "1: Insert key and value" << endl;
	cout << "2: Delete key" << endl;
	cout << "3: Get value[key]" << endl;
	cout << "4: Print" << endl;
	cout << "5: Delete BST" << endl;
	cout << "6: Get min. key" << endl;
	cout << "7: Get rank of a key" << endl;
	cout << "0: Exit" << endl;
	int option = 99;
	while (cin >> option)
	{
		if (0 == option)
			break;
		else if (1 == option)
		{
			int key;
			string value;
			cin >> key >> value;
			InsertKey(key, value);
		}
		else if (2 == option)
		{
			int key;
			cin >> key;
			DeleteKey(key);
		}
		else if (3 == option)
		{
			int key;
			cin >> key;
			cout << "Get: Key: " << key << " Value: " << GetValue(key) << endl;
		}
		else if (4 == option)
			PrintBst();
		else if (5 == option)
			DeleteBst();
		else if (6 == option)
			cout << "Get Min: Key: " << GetMin() << endl;
		else if (7 == option)
		{
			int key;
			cin >> key;
			cout << "Get: Rank: " << GetRank(key) << " Key: " << key << endl;
		}
		else
		{
			cout << "Invalid option: " << option << endl;
			break;
		}
	}
	DeleteBst();
	*/

	//Graph
	const int N = 15;
	int TestConnect[N][2] = { { 0,1 }, { 0,2 }, { 3,4 }, { 5,8 }, { 4,0 },
							  { 8,8 }, { 9,2 }, { 6,4 }, { 2,5 }, { 0,5 },
							  { 3,0 }, { 7,8 }, { 2,7 }, { 9,4 }, { 6,8 } };

	for (int i = 0; i < N; ++i)
	{
		MyDfs.AddConnection(TestConnect[i][0], TestConnect[i][1]);
		MyBfs.AddConnection(TestConnect[i][0], TestConnect[i][1]);
	}
	MyDfs.PrintConnections();
	MyDfs.PerformDfs();

	MyBfs.PrintConnections();
	MyBfs.PerformBfs();
    
    return 0;
}

void InsertKey(const int key, const string value)
{
	Head = InsertKey(Head, key, value);
	cout << "Insert complete: Key: " << key << " Value: " << value << endl;
}

shared_ptr<BstType> InsertKey(shared_ptr<BstType> Node, const int key, const string value)
{
	//Key is not present. Create a new node, and update size
	if (nullptr == Node)
	{
		shared_ptr<BstType> NewNode = make_shared<BstType>();
		NewNode->key = key;
		NewNode->value = value;
		return NewNode;
	}

	if (key == Node->key)
		//Key is already present. Just replace value with new one and return;
		Node->value = value;
	else if (key > Node->key)
		Node->right = InsertKey(Node->right, key, value);
	else
		Node->left = InsertKey(Node->left, key, value);
	
	Node->size = GetSize(Node->left) + GetSize(Node->right) + 1;
	return Node;
}

void DeleteKey(const int key)
{
	Head = DeleteKey(Head, key);
}

shared_ptr<BstType> DeleteKey(shared_ptr<BstType> Node, const int key)
{
	if (nullptr == Node)
	{
		cout << "Delete: Key: " << key << " not found" << endl;
		return Node;
	}

	if (key > Node->key)
		Node->right = DeleteKey(Node->right, key);
	else if (key < Node->key)
		Node->left = DeleteKey(Node->left, key);
	else
	{
		cout << "Delete: Key: " << key << endl;

		//Key found
		if (nullptr == Node->left)
		{
			shared_ptr<BstType> ret = Node->right;
			Node.reset();
			return ret;
		}
		else if (nullptr == Node->right)
		{
			shared_ptr<BstType> ret = Node->left;
			Node.reset();
			return ret;
		}
		else
		{
			//Since key to be deleted has 2 childs, Replace the node to be deleted with its min. node of right tree
			shared_ptr<BstType> MinNode = GetMin(Node->right);
			Node->key = MinNode->key;
			Node->value = MinNode->value;
			Node->right = DeleteMin(Node->right); //Delete the min. node of right tree
		}
	}

	Node->size = GetSize(Node->left) + GetSize(Node->right) + 1;
	return Node;
}

string GetValue(const int key)
{
	shared_ptr<BstType> temp = Head;
	while (nullptr!=temp)
	{
		if (key == temp->key)
			return temp->value;
		else if (key > temp->key)
			temp = temp->right;
		else
			temp = temp->left;
	}
	return "Not Found!!!";
}

void DeleteBst(void)
{
	DeleteBst(Head);
	Head = nullptr;
	cout << "Delete BST complete" << endl;
}

void DeleteBst(shared_ptr<BstType> Node)
{
	if (nullptr == Node)
		return;
	DeleteBst(Node->left);
	DeleteBst(Node->right);
	Node.reset();
}

void PrintBst(void)
{
	PrintBst(Head);
	cout << "Print BST complete" << endl;
}

void PrintBst(const shared_ptr<BstType> Node)
{
	if (nullptr == Node)
		return;
	PrintBst(Node->left);
	cout << "Print: Key: " << Node->key << " Value: " << Node->value << " Size: " << Node->size << endl;
	PrintBst(Node->right);
}

int GetSize(const shared_ptr<BstType> Node)
{
	if (nullptr == Node)
		return 0;
	return Node->size;
}

int GetRank(const int key)
{
	return GetRank(Head, key);
}

int GetRank(const shared_ptr<BstType> Node, const int key)
{
	if (nullptr == Node)
		return 0;

	if (key == Node->key)
		return GetSize(Node->left);
	else if (key < Node->key)
		return GetRank(Node->left, key);
	else
	{
		return (GetRank(Node->right, key) + GetSize(Node->left)+1);
	}
}

int GetMin(void)
{
	if (nullptr != Head)
		return GetMin(Head)->key;
	else
		return -9999;
}

shared_ptr<BstType> GetMin(const shared_ptr<BstType> Node)
{
	if (nullptr == Node->left)
		return Node;
	return GetMin(Node->left);
}

void DeleteMin(void)
{
	if(nullptr != Head)
		Head = DeleteMin(Head);
}

shared_ptr<BstType> DeleteMin(shared_ptr<BstType> Node)
{
	if (nullptr == Node->left)
	{
		shared_ptr<BstType> ret = Node->right;
		Node.reset();
		return ret;
	}		
	Node->left = DeleteMin(Node->left);
	Node->size = GetSize(Node->left) + GetSize(Node->right) + 1;
	return Node;
}

int GetMax(void)
{
	return 0;
}

int NumWaysNsum(const int Sum, const vector<int>::size_type VecIdx, unordered_map<string,int> &ResultMap, const vector<int> &InputVec)
{
    static int FuncCnt = 0;
  
    ++FuncCnt;
    cout << "Total: " << Sum << " Inputs: ";
    for(auto i = VecIdx; i < InputVec.size(); --i)
        cout << InputVec[i] << " ";
    cout << "FuncCnt: " << FuncCnt << endl;
    
    if(0 == Sum)
        return 1;
    else if (Sum < 0)
        return 0;
    else if(VecIdx >= InputVec.size())
        return 0;
    
    int Ret = 0;
    string key = to_string(Sum) + ':' + to_string(VecIdx);
    auto MapItr = ResultMap.find(key);
    if(MapItr != ResultMap.end())
        return MapItr->second;
    
    if(Sum < InputVec[VecIdx])
        Ret = NumWaysNsum(Sum, VecIdx-1, ResultMap, InputVec);
    else
        Ret = NumWaysNsum(Sum-InputVec[VecIdx], VecIdx-1, ResultMap, InputVec) + NumWaysNsum(Sum, VecIdx-1, ResultMap, InputVec);
    ResultMap[key] = Ret; //Create new element in Map and store the result
    
    cout << "Ret: Total: " << Sum << " Inputs: ";
    for(auto i = VecIdx; i < InputVec.size(); --i)
        cout << InputVec[i] << " ";
    cout << "NumWays: " << Ret << endl;
    
    return Ret;
}

int NumWaysStair(const int Total, unordered_map<int,int> &ResultMap, const vector<int> &Steps)
{
    static int FuncCnt = 0;
   
    if(0 == Total)
        return 1;
 
    auto MapItr = ResultMap.find(Total);
    if(MapItr != ResultMap.end())
        return MapItr->second;
    
    ++FuncCnt;
    cout << "Total: " << Total << " FuncCnt: " << FuncCnt << endl;
    
    int Ret = 0;
    for(auto &i:Steps)
        if(Total-i>=0)
            Ret += NumWaysStair(Total-i, ResultMap, Steps);
    
    ResultMap[Total] = Ret; //Create new element in Map and store the result
    
    cout << "Ret: Total: " << Total << " NumWays: " << Ret << endl;
    
    return Ret;
}

/*
void SetRoot(const int R, const int C, const int RootVal)
{
    if((R<0) || (R>=M)||
       (C<0) || (C>=N))
        return;
    if((RootVal == ComponentRoot[R][C]) || (!ComponentInputs[R][C]))
        return;
    
    //Set RootVal for current node and for all its 4 ngbrs
    ComponentRoot[R][C] = RootVal;
    for(auto &k : MatrixNgbrs)
        SetRoot(R+k.first, C+k.second, RootVal);
}
*/
    
void InsertHeapTree(vector<int> &heaptree, int val)
{
    //Add new value to the end of the tree, and then swim up to maintain heap
    heaptree.push_back(val);
    HeapSwimUp(heaptree, heaptree.size()-1);
}

int GetHeapTree(vector<int> &heaptree)
{
    //return 0th index which will be the max. After returning the 0th index,
    //copy the last element into 0th index, reduce the tree size and then
    //swim down to maintain heap
    if(!heaptree.size())
        return -9999;
    
    int ret = heaptree[0];
    SwapIvec(heaptree,0,heaptree.size()-1);
    heaptree.pop_back();
    if(heaptree.size())
        HeapSwimDown(heaptree,0,heaptree.size());
    
    return ret;
}

void HeapSwimUp(vector<int> &heaptree, vector<int>::size_type k)
{
    //Parent of k will be (k-1)/2
    //Keep swapping until parent > child
    while((k>0) && (k<heaptree.size()))
    {
        if(heaptree[(k-1)/2] < heaptree[k])
        {
            SwapIvec(heaptree,k,(k-1)/2);
            k = (k-1)/2;
        }
        else
            break;
    }
}

void HeapSwimDown(vector<int> &heaptree, vector<int>::size_type k, vector<int>::size_type vecsize)
{
    //Children of k will be (2k+1) and (2K+2)
    //Keep swapping until parent > max(child1,child2)
    while(2*k+1 < vecsize)
    {
        auto c1 = 2*k+1, max = c1;
        auto c2 = (c1+1 < vecsize) ? c1+1 : c1;
        if(heaptree[c1] < heaptree[c2])
            max = c2;
        if(heaptree[max] > heaptree[k])
        {
            SwapIvec(heaptree,k,max);
            k = max;
        }
        else
            break;
    }
}

void ShuffleVec(vector<int> &ivec)
{
    uniform_int_distribution<unsigned> u(0,ivec.size()-1);
    //default_random_engine e_for_seed; //Use a random generator to generate a random seed
    //default_random_engine e_for_index(e_for_seed()); //Use the random seed to generate another random generator for randomizing
                                            //array index
    default_random_engine e_for_index;
    
    for(vector<int>::size_type i = 0; i < ivec.size(); ++i)
    {
        //Swap index i with random index u(e)
        int randindex1 = u(e_for_index);
        int randindex2 = u(e_for_index);
        //Sanity chk
        if((randindex1 < ivec.size()) && (randindex2 < ivec.size()))
            SwapIvec(ivec,randindex1,randindex2);
        //PrintSortVec(ivec);
    }
}

void MyQuickSort(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    static int reccnt = 0;
    cout << "Sort lo: " << lo << " hi: " << hi << " cnt: " << ++reccnt << endl;
    if((lo>=hi) || (lo>=ivec.size()) || (hi>=ivec.size()))
        return;
    auto j = MyQuickPartition(ivec,lo,hi);
    MyQuickSort(ivec,lo,j-1);
    MyQuickSort(ivec,j+1,hi);
    
    PrintSortVec(ivec);
}

vector<int>::size_type MyQuickPartition(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    static int reccnt = 0;
    cout << "                          Partition lo: " << lo << " hi: " << hi << " cnt: " << ++reccnt << endl;
    auto i = lo, j = hi+1;
    int v = ivec[lo]; //partition element
    
    while(true)
    {
        while((v > ivec[++i]) && (i <= hi)) {}
        while((v < ivec[--j]) && (j > lo)) {}
        if(i>=j)
            break;
        SwapIvec(ivec,i,j);
    }
    SwapIvec(ivec,lo,j);
    
    return j;
}

void SwapIvec(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    if(lo==hi)
        return;
    int temp = ivec[lo];
    ivec[lo] = ivec[hi];
    ivec[hi] = temp;
}

void MyMergeSort(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    static int reccnt = 0;
    cout << "Sort lo: " << lo << " hi: " << hi << " cnt: " << ++reccnt << endl;
    if(lo>=hi)
        return;
    auto mid = (lo+hi)/2;
    MyMergeSort(ivec,lo,mid);
    MyMergeSort(ivec,mid+1,hi);
    MyMerge(ivec,lo,hi);
    
    PrintSortVec(ivec);
}

void MyMerge(vector<int> &ivec, vector<int>::size_type lo, vector<int>::size_type hi)
{
    static int reccnt = 0;
    cout << "                          Merge lo: " << lo << " hi: " << hi << " cnt: " << ++reccnt << endl;
    
    if(ivec.size()<2)
        return;
    auto mid = (lo+hi)/2;
    
    vector<int> tempvec(ivec);
    
    decltype(lo) i = lo; //loop ctr for 1st half of the array. i.e lo to mid
    decltype(lo) j = mid+1; //loop ctr for 2nd half of the array. i.e mid+1 to hi
    
    for(decltype(lo) k = lo; k <= hi; ++k)
        if(j>hi)
            ivec[k] = tempvec[i++];
        else if(i>mid)
            ivec[k] = tempvec[j++];
        else if(tempvec[i] < tempvec[j])
            ivec[k] = tempvec[i++];
        else
            ivec[k] = tempvec[j++];
}

/*
bool operator< (const string &s1, const string &s2)
{
    if(s1.size()<s2.size())
        return true;
    else if(s1.size()>s2.size())
        return false;
    else
    {
        for(int i = 0; i < s1.size(); ++i)
            if (s1[i] < s2[i])
                return true;
            else if (s1[i] > s2[i])
                return false;
        return false;
    }
}
 */

void PrintSetlist(const list<unordered_set<int> > &setlist)
{
    int setcnt = 0;
    for(const auto &i : setlist)
    {
        ++setcnt;
        cout << "Set#" << setcnt << " ";
        for(const auto &j : static_cast<unordered_set<int>>(i))
        {
            cout << j << " ";
        }
        cout << endl;
    }
}

void PrintSortVec(const vector<int> &ivec)
{
    for(const auto &i : ivec)
        cout << i << " ";
    cout << endl;
}

template <typename T>
bool Popstack(stack<T> &S, T &Val)
{
    if(!S.empty())
    {
        Val = S.top();
        S.pop();
        return true;
    }
    return false;
}

long long int Factorial(const int N)
{
    static int cnt = 0;
    ++cnt;
    
    cout << "N: " << N << " Cnt: " << cnt << endl;
    
    if(N<2)
        return 1;
    else
        return (N*Factorial(N-1));
}

int BinSearchRecursive(const int S, const vector<int> &A, const vector<int>::size_type lo, const vector<int>::size_type hi)
{
    static int cnt = 0;
    ++cnt;
    
    decltype(lo) mid = (lo+hi)/2;
    if(lo>hi)
    {
        cout << "Case 1: " << lo << " " << mid << " " << hi << " " << cnt << endl;
        return -1;
    }
    if(A[mid] == S)
    {
        cout << "Case 2: " << lo << " " << mid << " " << hi << " " << cnt << endl;
        return static_cast<int>(mid);
    }
    else if(S < A[mid])
    {
        cout << "Case 3: " << lo << " " << mid << " " << hi << " " << cnt << endl;
        return BinSearchRecursive(S,A,lo,mid-1);
    }
    else
    {
        cout << "Case 4: " << lo << " " << mid << " " << hi << " " << cnt << endl;
        return BinSearchRecursive(S,A,mid+1,hi);
    }
}

void WriteToQueue(queue<string> &WordQ, const string &Word)
{
    if( WordQ.size() >= 6 )
        WordQ.pop();
    WordQ.push(Word);
}

int f(void)
{
    cout << "f" << endl;
    return 1;
}
int g(void)
{
    cout << "g" << endl;
    return 1;
}
int h(void)
{
    cout << "h" << endl;
    return 1;
}
int j(void)
{
    cout << "j" << endl;
    return 1;
}

bool BinSearch(const int i, vector<int> *pvec)
{
    auto low = pvec->begin();
    auto hi = pvec->end();
    auto mid = pvec->begin() + pvec->size()/2;
    
    while( low != hi )
    {
        if(i == *mid)
            return true;
        else if (i < *mid)
            hi = mid;
        else
            low = mid + 1;
        mid = low + (hi - low)/2;
    }
    return false;
}
