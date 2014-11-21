// worm.cpp : 定义控制台应用程序的入口点。
//
#include "stdafx.h"
using namespace cv; 
using namespace std;

const int COLOR_RANGE = 256;
const int BLANK = 255;
const int BLACK = 0;
const double INF = 1e100;
const double MIN_WORM_SIZE_RATE = 0.1;
const double MIN_WORM_SEG_RATE = 0.15;

typedef vector<int> vi;
typedef vector<vi> vii;

typedef vector<uchar> vu;
typedef vector<vu> vuu;

inline double sqr(double x)
{
	return x * x;
}

const int totalDir = 8;
const int dir[totalDir][2] = {{-1, -1},{-1, 0}, {-1,1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}};

template <class T>
inline void initVec(vector <vector <T> > &x, int n, int m, T v)
{
	x.resize(n);
	for (int i = 0; i < n; i++) x[i].resize(m, v);
}

class Worm
{
public:
	Worm(int s): size(s)
	{
	}
	vector <Point> margin;
	int size;
};

class Edge
{
public:
	Edge(int x, int y, int length) :x(x), y(y), length(length), alive(true), father(-1)
	{
	}
	int length, x, y;
	int father;
	bool alive;
};

typedef set <int> Edgelist;

class Graph
{
public:
	int n,m;
	vector <Edgelist> eList;
	vector <Edge> edges;
	vector <Point> nodes;
	int addNode(int x, int y)
	{
		nodes.push_back(Point(x, y));
		n++;
		eList.push_back(Edgelist());
		return nodes.size() - 1;
	}
	int addEdge(int x, int y, int length)
	{
		edges.push_back(Edge(x,y,length));
		int num = edges.size() - 1;
		eList[x].insert(num);
		eList[y].insert(num);
		m++;
		return num;
	}
	int addEdge(Edge e)
	{
		edges.push_back(e);
		int num = edges.size() - 1;
		eList[e.x].insert(num);
		eList[e.y].insert(num);
		m++;
		return num;
	}
	void removeEdge(Edge &tedge)
	{
		if (tedge.alive)
		{
			tedge.alive = false;
			m--;
		}
	}
	void unionNode(int x, int y)
	{
		printf("%d %d\n", x, y);
		eList[x].insert(eList[y].begin(), eList[y].end());
		for (vector<Edge>::iterator it = edges.begin(); it != edges.end(); it++)
		{
			if (it -> x == y)
			{
				it -> x = x;
				if (it -> y == x) removeEdge(*it);
			}
			if (it -> y == y)
			{
				it -> y = x;
				if (it -> x == x) removeEdge(*it);
			}
		}
		n--;
	}
	Graph(): n(0), m(0)
	{
		edges.clear();
		eList.clear();
		nodes.clear();
	}
	void clear()
	{
		n = 0;
		m = 0;
		eList.clear();
		edges.clear();
		nodes.clear();
	}
	void cleanEdge(int lowest)
	{
		for (vector <Edge>::iterator it = edges.begin(); it != edges.end(); it++)
		{
			if ( it -> alive && it-> length < lowest)
			{

				unionNode(it->x, it -> y);
				removeEdge(*it);
			}
		}
	
	}
	void printInfo()
	{
		printf("n = %d m = %d", n, m);
		printf("nodes:\n");
		for (vector <Point>::iterator it = nodes.begin(); it != nodes.end(); it++)
			printf("x:%d ,y:%d\n", it -> x, it -> y);
		printf("edges:\n");
		for (vector <Edge>::iterator it = edges.begin(); it != edges.end(); it++)
			if (it -> alive)
				printf("x:%d ,y:%d\n", it -> x, it -> y);
		printf("List:\n");
		for (int i = 0; i < eList.size(); i++)
		{
			printf("node %d:", i);
			for (Edgelist::iterator it = eList[i].begin(); it != eList[i].end(); it++)
				if (edges[*it].alive)
					printf("%d ", *it);
			printf("\n");
		}
	}
private:

};

class WormMat :public Mat
{

public:
	vector <int> region;
	vector <Worm> worms;
	vector< vector <uchar>> sign;
	int countNumber()
	{
		Graph g;
		vii visit;
		initVec(visit, rows, cols, 0);
		int ret = 0;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				if (!visit[i][j] && this -> at<uchar>(i,j) == BLACK)
				{
					g.clear();
					int maxLength = 0;

					constructGraph(*this, g, visit, i, j, Edge(-1,-1,0), maxLength, false);
					g.cleanEdge(maxLength * MIN_WORM_SEG_RATE);
					g.printInfo();
					ret += countNumber(g);
				}
		return ret;
	}

	WormMat () : Mat()
	{
		regionIsSet = 0;
		this->sign.clear();
	}
	WormMat (const Mat &b) : Mat(b)
	{
		regionIsSet = 0;
		sign.clear();
	}
	~WormMat()
	{
		~Mat();
	}
	WormMat extract(int layer, int s, int e)
	{
		s--;
		e--;
		if (regionIsSet != layer) countRegion(layer);
		WormMat ret = *this;
		for (Mat_<uchar>::iterator it = ret.begin<uchar>(); it!= ret.end<uchar>(); ++it)
		{
			int t = 0;
			while (*it > region[t + 1]) t++;
			if (t < s || t > e) *it = 255;
		}
		return ret;
	}

	void countRegion(int layer)
	{
		regionIsSet = layer;
		int color[COLOR_RANGE], tot[COLOR_RANGE];
		vii pre;
		vector <vector <double> >opt;
		region.resize(layer + 1, 0);
		initVec(pre, COLOR_RANGE, layer + 1, -1);
		initVec(opt, COLOR_RANGE, layer + 1, 0.0);
		double sumx[COLOR_RANGE], sumx2[COLOR_RANGE];
		memset(color, 0, sizeof color);

		for (Mat_<uchar>::iterator it= this -> begin<uchar>(); it!= this -> end<uchar>(); ++it)
			color[*it]++;
		/*for (int i = 0; i < COLOR_RANGE; i++)
		{
			printf("%d %d\n", i, color[i]);
		}*/
		pretreatment(tot, sumx, sumx2, color);
		for (int i = 1; i < COLOR_RANGE; i++)
			opt[i][1] = tot[i] ? (sumx2[i] / tot[i] - sqr(sumx[i] / tot[i])) : 0;
		for (int i = 1; i < COLOR_RANGE; i++)
			for (int k = 2; k <= layer; k++) opt[i][k] = INF;
		for (int i = 2; i < COLOR_RANGE; i++)
			for (int j = 1; j < i; j++)
			{
				double variance = 0;
				int n = tot[i] - tot[j - 1];
				if (n > 0)
				{
					double ex2 = (sumx2[i] - sumx2[j - 1]) / n, e2x = sqr((sumx[i] - sumx[j - 1]) / n);
					variance = ex2 - e2x;
				}
				for (int k = 2; k <= layer; k++)
				{
					if (opt[j][k - 1] + variance < opt[i][k])
					{
						opt[i][k] = opt[j][k - 1] + variance;
						pre[i][k] = j;
					}
				}
			}
			//printf("%d\n",pre[COLOR_RANGE - 1][layer]);
			for (int p = COLOR_RANGE - 1, i = layer; i > 1; i--, p = pre[p][i])
				region[i - 1] = pre[p][i];
			region[0] = 0;
			region[layer] = COLOR_RANGE - 1;
			for (int i = 0; i <= layer; i++) printf("%d ",region[i]);
	}

	void thinning()
	{
		static const int erasetable[256]={
0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,1,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,0,1,1,0,0,1,0,0,0
		};
		bool Finished = false;
		while(!Finished)
		{
			Finished = true;
			for (int y = 1; y < cols - 1; y++)
			{ 
				int x = 1; 
				while(x< rows - 1)
				{
					if(this -> at<uchar>(x,y) == BLACK)
					{
						if( (this -> at<uchar>(x - 1,y)==BLANK)|| (this -> at<uchar>(x + 1,y)==BLANK))
						{
							int num = 0;
							for (int i = 0; i < totalDir; i++)
								num |= ((this -> at<uchar>(x + dir[i][0], y + dir[i][1])) / 255) << i;
							if(erasetable[num] == 1)
							{
								this -> at<uchar>(x,y) = BLANK;
								Finished = false;
								x++;
							}
						}
					}
					x++;
				}
			}
			for (int x = 1;x< rows - 1; x++)
			{ 
				int y = 1;
				while(y < cols - 1)
				{
					if(this -> at<uchar>(x,y)== BLACK)
					{
						if( (this -> at<uchar>(x,y - 1) == BLANK)|| (this -> at<uchar>(x,y + 1) == BLANK))
						{
							int num = 0;
							for (int i = 0; i < totalDir; i++)
								num |= ((this -> at<uchar>(x + dir[i][0], y + dir[i][1])) / 255) << i;
							if(erasetable[num] == 1)
							{
								this -> at<uchar>(x,y) = BLANK;
								Finished = false;
								y++;
							}
						}
					}
					y++;
				}
			} 
		}
	}

	void split(WormMat &afterClean, WormMat &margin)
	{
		Mat &a = *this;
		int n = a.rows, m = a.cols;
		initVec(sign ,n,m, uchar(0));
		int cnt = 0, mSize = 0, tar = 0; 
		vector <int> size;
		size.clear();
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				if (!sign[i][j] && a.at<uchar>(i,j) < BLANK)
				{
					int nowSize = 0;
					myFloodFill(a, i, j, ++cnt, sign, nowSize);
					size.push_back(nowSize);
					if (nowSize > mSize)
					{
						mSize = nowSize;
						tar = cnt;
					}
				}
		cnt = 0;
		vii isMargin;
		initVec(isMargin, rows, cols, 0);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				if (sign[i][j] > cnt)
				{
					//printf("%d %d %d\n",size[sign[i][j] - 1], mSize, sign[i][j]);
					if (size[sign[i][j] - 1] < mSize * MIN_WORM_SIZE_RATE) fillInt(sign, i, j, 0);
					else
					{
						size[++cnt - 1] = size[sign[i][j] - 1];
						fillInt(sign, i, j, cnt);
						findMargin(sign, i, j, isMargin, 0);
					}
				}
		afterClean = Mat(rows, cols, CV_8UC1);
		margin = Mat(rows, cols, CV_8UC1);
		Mat_<uchar>::iterator it= afterClean.begin<uchar>();
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				if (sign[i][j]) *it++ = BLACK;
				else *it++ = BLANK;
		it= margin.begin<uchar>();
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				if (isMargin[i][j])
					*it++ = BLACK; 
				else *it++ = BLANK;
	}

private:
	int regionIsSet;

	void constructGraph(Mat &a, Graph &g, vii &visit, int x, int y, Edge e, int &maxLength, bool newNode)
	{
		e.length++;
		maxLength = max(maxLength, e.length);
		int cnt = 0;
		for (int i = 0; i < totalDir; i++)
		{
			int tx = x + dir[i][0];
			int ty = y + dir[i][1];
			if (tx >= 0 && tx <rows && ty >= 0 && ty < cols && a.at<uchar>(tx,ty) == BLACK)
			{
				if (!newNode && visit[tx][ty] > 1)
				{
					e.y = visit[tx][ty] - 2;
					g.addEdge(e);
					return;
					//return;
				}
				else
					if (!visit[tx][ty]) cnt++;
			}
		}
		newNode = cnt == 0 || cnt > 1 || e.x == -1;
		if (newNode)
		{
			e.y = g.addNode(x,y);
			//e.length = dist(g)
			if (e.x != -1) g.addEdge(e);
			visit[x][y] = e.y + 2;
			if (cnt == 0) return;
			e.x = e.y;
			e.y = -1;
			e.length = 0;
		}
		else
			visit[x][y] = 1;
		for (int i = 0; i < totalDir; i++)
		{
			int tx = x + dir[i][0];
			int ty = y + dir[i][1];
			if (tx >= 0 && tx <rows && ty >= 0 && ty < cols && a.at<uchar>(tx,ty) == BLACK && !visit[tx][ty]) constructGraph(a, g, visit, tx, ty, e, maxLength, newNode);
		}
	}

	int countNumber(Graph &g)
	{

		return 0;
	}

	void pretreatment(int *tot, double *sumx, double *sumx2, int *color)
	{
		memset(tot, 0, COLOR_RANGE * sizeof(int));
		memset(sumx, 0, COLOR_RANGE * sizeof(int));
		memset(sumx2, 0, COLOR_RANGE * sizeof(int));
		tot[0] = color[0];
		sumx[0] = 0;
		sumx2[0] = 0;
		for (int i = 1; i < COLOR_RANGE; i++)
		{
			tot[i] = tot[i - 1] + color[i];
			sumx[i] = sumx[i - 1] + double(color[i]) * i;
			sumx2[i] = sumx2[i - 1] + double(color[i]) * i * i;
		}
	}

	void myFloodFill(Mat &a, int x, int y, int cnt, vuu &sign, int &size)
	{
		sign[x][y] = cnt;
		//printf("%d %d\n", a.at<uchar>(x,y));
		size++;
		for (int i = 0; i < totalDir; i++)
		{
			int tx = x + dir[i][0];
			int ty = y + dir[i][1];
			if (tx >= 0 && tx <a.rows && ty >= 0 && ty < a.cols && sign[tx][ty] == 0 && a.at<uchar>(tx,ty) < BLANK)
				myFloodFill(a, tx, ty, cnt, sign, size);
		}
	}

	void fillInt(vuu &sign, int x, int y, int cnt)
	{
		sign[x][y] = cnt;
		for (int i = 0; i < totalDir; i++)
		{
			int tx = x + dir[i][0];
			int ty = y + dir[i][1];
			if (tx >= 0 && tx <rows && ty >= 0 && ty < cols && sign[tx][ty] != cnt && sign[tx][ty]) fillInt(sign, tx, ty, cnt);
		}
	}

	void findMargin(vuu &sign, int x, int y, vii &isMargin, int nowdir)
	{
		//printf("%d %d %d\n", x, y, isMargin[x][y]);
		for (int i = nowdir, j = 0; j < totalDir; j++)
		{
			int tx = x + dir[i][0];
			int ty = y + dir[i][1];
			if (tx >= 0 && tx < rows && ty >= 0 && ty < cols && sign[tx][ty] == sign[x][y])
			{
				if (!(isMargin[x][y] & (1 << i)))
				{
					isMargin[x][y] |= 1 << i;
					findMargin(sign, tx, ty, isMargin, (i + totalDir - 1) % totalDir);
				}
				break;
			}
			i = (i + 1) % totalDir;
		}
	}
};

Mat orgImg, imgGray;

void init(String fn)
{
	orgImg = imread(fn); 
	if (orgImg.empty()) 
	{
		fprintf(stderr,"Error: load image failed."); 
		exit(-1); 
	} 
	cvtColor(orgImg, imgGray, CV_RGB2GRAY);
	printf("%d\n", imgGray.type());

}

int _tmain(int argc, _TCHAR* argv[])
{
	init("1_3.png");
	WormMat	img = imgGray;

	imwrite("gray.bmp",img);
	medianBlur(img, img, 9);
	imwrite("blur.bmp",img);
	img.countRegion(5);
	WormMat newimg = img.extract(5,1,3);
	namedWindow("image", CV_WINDOW_AUTOSIZE);
	namedWindow("image2", CV_WINDOW_AUTOSIZE);

	namedWindow("image3", CV_WINDOW_AUTOSIZE);
	imshow("image2", newimg);
	imwrite("extract.bmp",newimg);
	WormMat afterClean, margin;
	newimg.split(afterClean, margin);
	imwrite("clean.bmp",afterClean);

	imshow("image", afterClean);
	afterClean.thinning();

	imwrite("thinning.bmp",afterClean);
	afterClean.countNumber();
	imshow("image3", afterClean);
	//imwrite("test2.bmp",afterClean);
	imshow("image2", margin);
	waitKey(); 
	return 0;
}

