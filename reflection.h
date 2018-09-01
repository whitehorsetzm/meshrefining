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
class Refletion
{
public:
    Refletion();
    Refletion(Refletion *a){
        this->subject_table=a->subject_table;
        this->curves=a->curves;
        this->sufaces=a->sufaces;
    }
    void initial(char *gm3file,HYBRID_MESH& mesh);//表的初始化可以通过subject()函数获取一个点离每个几何面的距离（或者点离一个几何面的环的距离），选取最近的。最后我们要通过点与面的关系得到
    //面与面的关系，于是我们需要一个数据结构来保存这个结果等待处理。
    Vector subject_test(int patch_ID_1,int patch_ID_2,Vector coord);
    void attach_face(int face_ID,int patch_ID);
    int detach_face(int face_ID);
    Vector subject_face_id(int face_ID_1,int face_ID_2,Vector coord);
    Vector subject_patch_id(int patch_ID_1,int patch_ID_2,Vector coord);
    FergusonCurve* findcurve(int patch_id_1,int patch_id_2,Vector coord);
    int findsurface(int curve_id_1,int curve_id_2,int curve_id_3,Vector vector);
    void edge_to_face(HYBRID_MESH &mesh,int &face_id_1, int &face_id_2,string edge_name);
    void initial_edge_table(HYBRID_MESH &mesh);
     map<int,int> subject_table;
    FergusonSurface *sufaces;
    FergusonCurve   *curves;
    vector<string> lines;

};
//extern table _table;
class temp{
public:
    int patch_id;
    int curve_id;
    int node_id;
    temp():patch_id(-1),curve_id(-1),node_id(-1){/*cout<<"sucess"<<endl;*/}
};

#endif // TABLE_H
