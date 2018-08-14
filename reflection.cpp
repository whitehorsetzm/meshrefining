#include "reflection.h"
#include "float.h"
#include "fstream"
#include "strstream"
table _table;
FergusonSurface *data_surface;
FergusonCurve *data_curve;
GBSolid gbsolid;
DiscretSolid discretsolid;
HYBRID_MESH mesh;
temp *ref_table;
map<int,int> face_patch;
void datainitail(char *gm3file,HYBRID_MESH& mesh){
    read_gm3(gm3file,&gbsolid);
    vertex_G=gbsolid.vertex_G;
    GBCurve *curves=gbsolid.curve_G;
    GBFace *faces=gbsolid.face_G;
    int nf=gbsolid.nloops_G;
    int nc=gbsolid.ncurves_G;
    data_surface = new  FergusonSurface[nf];
    data_curve = new  FergusonCurve[nc];
    for(int i=0;i<nc;++i){
    createFergusonCurve(&gbsolid,&curves[i],&data_curve[i]);
    }
    for(int i=0;i<nf;++i){
    createFergusonFace(faces[i],&data_surface[i]);
    }
    double u,v;
    Vector coord;
    double distance,min,cls_u;
    double vol=0.0001;
    int k,f[nf],m=0;
    ref_table=new temp[mesh.NumNodes];





//    // suanfa 1
//    vertex *_vertex=new vertex[mesh.NumNodes];
//    for(int i=0;i<mesh.NumTris;++i){
//        for(int j=0;j<3;++j)
//        {
//    //        cout<<mesh.pTris[i].vertices[j]<<endl;
//            _vertex[mesh.pTris[i].vertices[j]].neighbor.insert(mesh.pTris[i].vertices[j+1%3]);
//            _vertex[mesh.pTris[i].vertices[j]].neighbor.insert(mesh.pTris[i].vertices[j+2%3]);
//        }
//    }

//    vector<vector<int>> curve_node;
//    for(int i=0;i<nc;++i){
//        vector<int> temp;
//        for(int j=0;j<mesh.NumNodes;++j){
//            data_curve[i].project(mesh.nodes[j].coord,&coord,&cls_u);
//            distance=mesh.nodes[j].coord.getDistance(coord);
//            if(distance<vol){
//                temp.push_back(j);
//                if(_vertex[j].on_curve==false){
//                    _vertex[j].curve=&data_curve[i];
//                    _vertex[j].on_curve=true;
//                    _vertex[j].flag=i;
//                }
//             }
//        }
//        curve_node.push_back(temp);
//   //     cout<<"i="<<i<<"   "<<temp.size()<<endl;
//    }

////    for(int i=0;i<curve_node.size();++i){
////        for(int j=0;j<curve_node[i].size();++j){
////                 cout<<"curve_ID= "<<i<<"node_ID= "<<curve_node[i][j]<<endl;
////        }
////    }
//    for(int i=0;i<curve_node.size();++i){
//    for(int j=0;j<curve_node[i].size();++j){
//    //      cout<<curve_node[i][j]<<"  ";
//    }
//   // cout<<endl;
//    }
//     vector<set<int>> loop_node;
//    for(int i=0;i<nf;++i){
//        set<int> temp;
//        for(int j=0;j<gbsolid.loop_G[i].nlc;++j){
//      //      cout<<"the "<<gbsolid.loop_G[i].cp_t[j]->index<<" curve "<<"has "<<curve_node[gbsolid.loop_G[i].cp_t[j]->index].size()<<endl;
//            for(int k=0;k<curve_node[gbsolid.loop_G[i].cp_t[j]->index].size();++k)
//            {
//                temp.insert(curve_node[gbsolid.loop_G[i].cp_t[j]->index][k]);
//              //    cout<<curve_node[gbsolid.loop_G[i].cp_t[j]->index][k]<<endl;
//            }

//        }
//        loop_node.push_back(temp);
//    }
//        for(int i=0;i<loop_node.size();++i){
//            for(auto j=loop_node[i].begin();j!=loop_node[i].end();++j){
//       //              cout<<"face_ID= "<<i<<"node_ID= "<<*j<<endl;
//            }
//        }


//      //  for(int i=0;i<mesh.NumNodes;++i)
//      //           cout<<"_vertex[i]"<<_vertex[i].neighbor.size()<<endl;
//for(int i=0;i<loop_node.size();++i){
//    for(auto j=loop_node[i].begin();j!=loop_node[i].end();++j){
//        for(auto k=_vertex[*j].neighbor.begin();k!=_vertex[*j].neighbor.end();++k)
//            {
//            data_surface[i].project(mesh.nodes[*k].coord,&u,&v);
//            data_surface[i].param_to_coord(u,v,&coord);
//            if(i=-52)
//                cout<<coord.getDistance(mesh.nodes[*k].coord)<<"?????????????????"<<endl;    //problem;
//            if(coord.getDistance(mesh.nodes[*k].coord)<vol&&_vertex[*k].flag==-1)
//            {
//                _vertex[*k].flag=i;
//                _vertex[*k].surface=&data_surface[i];
//            }
//        }
//   }
//}
//for(int i=0;i<mesh.NumNodes;++i){
//    cout<<_vertex[i].flag<<endl;
//}

//for(int i=0;i<mesh.NumNodes;++i){
//    if(_vertex[i].on_curve==true){
//         ref_table[i].curve_id=_vertex[i].flag;
//         ref_table[i].patch_id=-1;
//         ref_table[i].node_id=i;
//    }
//    else {
//        ref_table[i].patch_id=_vertex[i].flag;
//        ref_table[i].curve_id=-1;
//        ref_table[i].node_id=i;
//   }

//}





      //suanfa 2
    for(int i=0;i<mesh.NumNodes;++i){
        min==DBL_MAX;
        k=0;
        for(int j=0;j<nf;++j){
            data_surface[j].project(mesh.nodes[i].coord,&u,&v);
            data_surface[j].param_to_coord(u,v,&coord);
            distance=mesh.nodes[i].coord.getDistance(coord);
            if(distance<=vol){
                f[m]=j;
                m++;
                k=j;
            }
           if(distance<min){
                min=distance;
                k=j;
            }
        }
      //  cout<<"distance:"<<min<<endl;
      //   cout<<m<<endl; //distance has problem
        if(m>=2){
            min=DBL_MAX;
            for(int j=0;j<nc;++j){
                data_curve[j].project(mesh.nodes[i].coord,&coord,&cls_u);
                distance=mesh.nodes[i].coord.getDistance(coord);
                if(distance<min){
                     min=distance;
                     k=j;
                 }
            }
            data_curve[k].project(mesh.nodes[i].coord,&coord,&cls_u);
            distance=mesh.nodes[i].coord.getDistance(coord);
//             cout<<"k="<<k<<"curve:"<<distance<<endl;

            ref_table[i].curve_id=k;
            ref_table[i].node_id=i;
        }
        else{

            data_surface[k].project(mesh.nodes[i].coord,&u,&v);
            data_surface[k].param_to_coord(u,v,&coord);
            distance=mesh.nodes[i].coord.getDistance(coord);
          //  cout<<"k="<<k<<"face:"<<distance<<endl;

             ref_table[i].patch_id=k;
             ref_table[i].node_id=i;
        }
         m=0;
         min=DBL_MAX;
    }





//suan fa 3
vector<vector<int>> curve_head;
vector<int> temp;
for(int i=0;i<nc;++i){
    for(int j=0;j<mesh.NumNodes;++j){
        if(data_curve[i].start_coord().getDistance(mesh.nodes[j].coord)<0.00001){
            temp.push_back(j);
        }
        if(data_curve[i].end_coord().getDistance(mesh.nodes[j].coord)<0.00001){
            temp.push_back(j);
        }
    }
    curve_head.push_back(temp);
    temp.clear();
}
discretsolid.NumFacets=mesh.NumTris;
discretsolid.NumPoints=mesh.NumNodes;
discretsolid.discreteFacets=new DiscretFacet[discretsolid.NumFacets];
discretsolid.discretPoints=new DiscretPoint[discretsolid.NumPoints];
for(int i=0;i<discretsolid.NumFacets;++i){
    discretsolid.discretPoints[i].index=i;
           discretsolid.discretPoints[i].x=mesh.nodes[i].coord.x;
           discretsolid.discretPoints[i].y=mesh.nodes[i].coord.y;
           discretsolid.discretPoints[i].z=mesh.nodes[i].coord.z;
}


for(int i=0;i<discretsolid.NumFacets;++i){
    discretsolid.discreteFacets[i].index=i;
    discretsolid.discreteFacets[i].points[0]=mesh.pTris[i].vertices[0];
    discretsolid.discreteFacets[i].points[1]=mesh.pTris[i].vertices[1];
    discretsolid.discreteFacets[i].points[2]=mesh.pTris[i].vertices[2];
}
buildRelationshipByPoint(&discretsolid);
buildFacetRelationshipByEdge(&discretsolid);


for(int i=0;i<discretsolid.NumPoints;++i){
     if(discretsolid.discretPoints[i].linkedPoints.size()!=discretsolid.discretPoints[i].linkedFacets.size()){
         discretsolid.discretPoints[i].linkedPoints.erase(0);  //problem?
     }
}

vector<vector<int>> curve_node;
int start_node;
int end_node;
int pre_node=-1;
DiscretPoint current_node;
//for(auto j=discretsolid.discretPoints[0].linkedPoints.begin();j!=discretsolid.discretPoints[0].linkedPoints.end();++j){
//    cout<<*j<<endl;
//}
for(int i=0;i<nc;++i){
temp.clear();
start_node=curve_head[i][0];
end_node=curve_head[i][1];
//cout<<"start_node  "<<start_node<<"   end_node  "<<end_node<<endl;
current_node=discretsolid.discretPoints[start_node];
temp.push_back(start_node);
while(current_node.index!=end_node){
//cout<<current_node.index<<endl;
   for(auto j=current_node.linkedPoints.begin();j!=current_node.linkedPoints.end();++j){
      //               cout<<"neig  "<<*j<<endl;
       data_curve[i].project(mesh.nodes[*j].coord,&coord,&cls_u);
       distance=mesh.nodes[*j].coord.getDistance(coord);
    //    cout<<"pre "<<pre_node<<endl;
       if(distance<0.05&&(*j)!=pre_node){
       temp.push_back(*j);
  //     cout<<"add "<<*j<<endl;
       pre_node=current_node.index;
       current_node=discretsolid.discretPoints[*j];
       break;
   }
}
}
pre_node=-1;
curve_node.push_back(temp);
}

for(int i=0;i<curve_node.size();++i){
    cout<<"curve_id"<<i;
    for(int j=0;j<curve_node[i].size();++j)
    cout<<"node "<<curve_node[i][j]<<"  ";
    cout<<endl;
}
vector<set<int>> iner_loop;
vector<set<int>> loop;
set<int> _temp;
for(int i=0;i<gbsolid.nloops_G;++i){
    _temp.clear();
     for(int j=0;j<gbsolid.loop_G[i].nlc;++j){
        for(int k=0;k<curve_node[gbsolid.loop_G[i].cp_t[j]->index].size();++k){
            int n=curve_node[gbsolid.loop_G[i].cp_t[j]->index][k];
            _temp.insert(n);
        }
    }
    loop.push_back(_temp);
}

for(int i=0;i<gbsolid.nloops_G;++i){
    _temp.clear();
    for(auto j=loop[i].begin();j!=loop[i].end();++j){
        current_node=discretsolid.discretPoints[*j];
         for(auto k=current_node.linkedPoints.begin();k!=current_node.linkedPoints.end();++k)
          {
                if(loop[i].find(*k)==loop[i].end()){
             data_surface[i].project(mesh.nodes[*k].coord,&u,&v);
            data_surface[i].param_to_coord(u,v,&coord);
            distance=mesh.nodes[*k].coord.getDistance(coord);
            cout<<"loop="<<i<<"  node="<<*k<<"  distance= "<<distance<<endl;
            if(distance<0.05)
               {
                _temp.insert(*k);
               }
          }
         }
    }
    iner_loop.push_back(_temp);
}
for(int i=0;i<gbsolid.nloops_G;++i){
 cout<<"loop_size "<<i<<" "<<loop[i].size()<<endl;
 cout<<"iner_loop_size "<<i<<" "<<iner_loop[i].size()<<endl;
// for(auto j=loop[i].begin();j!=loop[i].end();++j){
//     cout<<" "<<*j;
// }
// cout<<endl;
}
for(int i=0;i<4;++i){
    cout<<"     "<<endl;
    for(auto j=iner_loop[i].begin();j!=iner_loop[i].end();++j){
        cout<<" "<<*j;
    }
}

    cout<<"reflection success!"<<endl;
    ofstream file;
    file.open("reflection_table.txt");
    if(file.is_open()){
        file<<"node_id"<<"  "<<"curve_id"<<"  "<<"patch_id"<<endl;
        for(int i=0;i<mesh.NumNodes;++i){
            file<<ref_table[i].node_id<<"  \t\t"<<ref_table[i].curve_id<<"  \t\t"<<ref_table[i].patch_id<<endl;
        }
        file.close();
    }

}  //#ok


table::table()
{

}
void table::initial(char *gm3file,HYBRID_MESH& mesh)
{
    datainitail(gm3file,mesh);
 //   ref_table=new temp[mesh.NumNodes];
//    while(file>>s)){
//        for(int i=0;i<3;++i)
//        {
//        }
//    }
    for(int i=0;i<mesh.NumTris;++i){
        int patch_id1=ref_table[mesh.pTris[i].vertices[0]].patch_id;
        int patch_id2=ref_table[mesh.pTris[i].vertices[1]].patch_id;
        int patch_id3=ref_table[mesh.pTris[i].vertices[2]].patch_id;
      //  cout<<patch_id1<<" "<<patch_id2<<" "<<patch_id3<<endl;
        if(patch_id1!=-1&&patch_id2!=-1&&patch_id3!=-1)
        {
            face_patch[i]=patch_id1;
        }
        else if(patch_id1==-1){
            if(patch_id2==-1)
                face_patch[i]=patch_id3;
            else if(patch_id3==-1)
                face_patch[i]=patch_id2;
            else
                face_patch[i]=patch_id2;
        }
        else if(patch_id2==-1){
            if(patch_id1==-1)
                face_patch[i]=patch_id3;
            else if(patch_id3==-1)
                face_patch[i]=patch_id1;
            else
            face_patch[i]=patch_id1;
        }
        else if(patch_id3==-1){
            if(patch_id1==-1)
                face_patch[i]=patch_id2;
            else if(patch_id2==-1)
                face_patch[i]=patch_id1;
            else
            face_patch[i]=patch_id1;
        }
    }

    for(int i=0;i<mesh.NumTris;++i){
        gbsolid.face_G[face_patch[i]].facets.insert(i);

    }
    this->curves=data_curve;
    this->sufaces=data_surface;
    set<int>::iterator it;
    for(int i=0;i<gbsolid.nloops_G;i++){
        it=gbsolid.face_G[i].facets.begin();
        for(int j=0;j<gbsolid.face_G[i].facets.size();++j){

            if(it!=gbsolid.face_G[i].facets.end())
           {

           // cout<<"patch_id "<<gbsolid.face_G[i].index<<"face id "<<*(it)<<endl;
            if(subject_table[*(it)]==0)
            subject_table[*(it)]=gbsolid.face_G[i].index+1;//problem
            it++;
           }
        }
        }

}

void table::attach_face(int face_ID, int patch_ID)
{
    subject_table[face_ID]=patch_ID;
}

int table::detach_face(int face_ID)
{
     int result =(*subject_table.find(face_ID)).second;
     subject_table.erase(face_ID);
     return result;
}
 int count=0;
Vector table::subject(int face_ID_1,int face_ID_2,Vector coord)
{
   Vector cls_pnt=coord;
   double u,v,cls_u;
   if(subject_table[face_ID_1]!=subject_table[face_ID_2])
   {
     //cout<<"face_ID_1: "<<face_ID_1<<"face_ID_1: "<<face_ID_2<<endl;
        FergusonCurve *temp;
      temp=findcurve(subject_table[face_ID_1]-1,subject_table[face_ID_2]-1,coord);
             if(temp==nullptr)
             {    cout<<"can't find curve!!!!!!!"<<endl;
                 cout<<"face_ID_1: "<<face_ID_1<<"face_ID_1: "<<face_ID_2<<endl;
                 cout<<"patch_id_1="<<subject_table[face_ID_1]<<"patch_id_2="<<subject_table[face_ID_2]<<endl;
             }
         else
              temp->project(coord,&cls_pnt,&cls_u);//problem

     //  cout<<"test";

   }
   else
   {
     //  cout<<"face_ID_1: "<<face_ID_1<<"face_ID_1: "<<face_ID_2<<endl;
     //  cout<<"patch_id_1="<<subject_table[face_ID_1]<<"patch_id_2="<<subject_table[face_ID_2]<<endl;
       this->sufaces[subject_table[face_ID_1]-1].project(coord,&u,&v);
       this->sufaces[subject_table[face_ID_1]-1].param_to_coord(u,v,&cls_pnt);
   }
   return cls_pnt;
}

FergusonCurve* table::findcurve(int patch_id_1, int patch_id_2,Vector coord)
{
    FergusonCurve *temp=new FergusonCurve;
  //  cout<<"patch_id_1="<<patch_id_1<<"patch_id_2="<<patch_id_2<<endl;
  //  cout<<gbsolid.loop_G[patch_id_1].nlc<<"  "<<gbsolid.loop_G[patch_id_2].nlc<<endl;
    for(int i=0;i<gbsolid.loop_G[patch_id_1].nlc;++i){
       for(int j=0;j<gbsolid.loop_G[patch_id_2].nlc;++j){
  //  cout<<"gbsolid.loop_G[patch_id_1].cp_t[i]="<<gbsolid.loop_G[patch_id_1].cp_t[i]->index<<"  gbsolid.loop_G[patch_id_2].cp_t[j]"<<gbsolid.loop_G[patch_id_2].cp_t[j]->index<<endl;
    if(gbsolid.loop_G[patch_id_1].cp_t[i]->index==gbsolid.loop_G[patch_id_2].cp_t[j]->index)
    {
        createFergusonCurve(&gbsolid,gbsolid.loop_G[patch_id_1].cp_t[i],temp);

        Vector cls_pnt;
        double cls_u;
        temp->project(coord,&cls_pnt,&cls_u);
        cout<<"distance = "<<cls_pnt.getDistance(coord)<<endl;
        if(cls_pnt.getDistance(coord)<0.05)         //problem;
               return temp;
    }
       }
    }
    return nullptr;
}
