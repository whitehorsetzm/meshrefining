#ifndef TABLE_H
#define TABLE_H
#include<map>
#include<vector>
#include<set>
#include<queue>
#include"dataclass.h"
#include"dataio.h"
#include"ferguson_curve.h"
#include"ferguson_surface.h"
#include"adaptive_io.h"
using namespace std;
using namespace EEMAS;
extern FergusonSurface *data_surface;
extern FergusonCurve *data_curve;
extern GBSolid gbsolid;
extern DiscretSolid dssolid;
extern HYBRID_MESH mesh;
void datainitail(char *gm3file,HYBRID_MESH& mesh);
class vertex{
public:
    vertex(){
        neighbor.clear();
        flag=-1;
        surface=nullptr;
        curve=nullptr;
        index=-1;
        face.clear();
        on_curve=false;
    }
  int index;
  set<int> neighbor;
  set<int> face;  //face contain this pointer
  FergusonSurface *surface;
  FergusonCurve   *curve;
  int flag;
  bool on_curve;
};
class table
{
public:
    table();
    table(table *a){
        this->subject_table=a->subject_table;
        this->curves=a->curves;
        this->sufaces=a->sufaces;
    }
    void initial(char *gm3file,HYBRID_MESH& mesh);//表的初始化可以通过subject()函数获取一个点离每个几何面的距离（或者点离一个几何面的环的距离），选取最近的。最后我们要通过点与面的关系得到
    //面与面的关系，于是我们需要一个数据结构来保存这个结果等待处理。
    void attach_face(int face_ID,int patch_ID);
    int detach_face(int face_ID);
    Vector subject(int face_ID_1,int face_ID_2,Vector coord);
    FergusonCurve* findcurve(int patch_id_1,int patch_id_2,Vector coord);
     map<int,int> subject_table;
    FergusonSurface *sufaces;
    FergusonCurve   *curves;
};
extern table _table;
class temp{
public:
    int patch_id;
    int curve_id;
    int node_id;
    temp():patch_id(-1),curve_id(-1),node_id(-1){/*cout<<"sucess"<<endl;*/}
};
//int vn=mesh.NumNodes;
//vertex* vertices=new vertex[vn];//下标为index
//void init(){
//    for(int i=0;i<elemsize;i++){
//        for(int j=0;j<3;++j){
//            vertices[mesh.pTris[i].vertices[j]].index=mesh.pTris[i].vertices[j];
//            vertices[mesh.pTris[i].vertices[j]].neighbor.insert(mesh.pTris[i].vertices[(j+1)%3]);
//            vertices[mesh.pTris[i].vertices[j]].neighbor.insert(mesh.pTris[i].vertices[(j+2)%3]);
//            vertices[mesh.pTris[i].vertices[j]].face.insert(mesh.pTris[i].index);
//        }
//    }
//}

//class loop_{
//  set<int>            node;//index
//  FergusonCurve     *curve;//loop由几段curve组成
//};
//class face_{
//  set<int>            node;
//  FergusonSurface    *suface;
//};
//loop_ edges[nloop];
//for(int i=0;i<n_node;++i){
//    node[i]
//}
//每个loop_里储存着几条几何曲线构成的环和这些曲线上所有的网格点.
//保存点与几何曲线的关系之后，可以得到网格的一圈顶点，根据这圈顶点作填色算法可以得到一个几何面对应的所有三角面。
//填色算法需要数据结构的支持
//该数据结构需要保存邻接点的信息，点的邻接信息可以通过遍历所有三角面得到。

//vertex vetices[vn] 下标为顶点的全局INDEX，遍历三角面片，对3个顶点进行操作，尝试向neighbor中加入另外两个点，若neighbor中没有这两个点，则直接加入，否则跳过。
//采用图算法遍历顶点时候，若当前顶点为边界点，则需要判断下一个顶点是否在需要的几何面上，如果是，则标记并访问，如果不是，则标记抛弃。
//保存边界点；
//如何计算某个点是否在某个几何面？ 将该点向几何面投影，然后将投影坐标与顶点坐标计算距离。（尝试）

//face_ face*=new face_[nf];
//for(int i=0;i<nloops;++i){
//    queue<int> q;
//    q.push(edges[i].node[0]);
//    while(!q.empty()){
//        int n=q.front();
//        face[i].node.insert(n);
//        q.pop();
//        vertices[n].flag=true;
//    for(int i=0;i<vertices[n].neighbor.size();++i){
//        if(!edges[i].node.find(t))
//        {
//            //如果找到，说明在边界上
//            //则要判断下个点是在面内还是面外
//            if(//面内&&vertices[n].flag==false)
//                    )
//               q.push(vertices[n].neighbor[i]);
//        else{
//            if(vertices[n].flag==false)
//            q.push(vertices[n].neighbor[i]);
//        }
//    }
//    }

//}
// }
//   face数组中保存了每个几何面所对应的网格点
//   edges数组中保存了每个几何线对应的网格点
//点全部处理完之后，要转化为面与面的关系。如何通过点找面？
//重新从mesh的所有面开始处理，此时我们有点属于哪个面或曲线的信息，这样处理起来就很快了。
//处理具体顺序，先对mesh所有面进行遍历，判断每个面中3个点是否有点属于某，


#endif // TABLE_H
