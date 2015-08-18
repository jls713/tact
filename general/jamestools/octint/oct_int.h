class node1{
	private:
		double loc;
		int level;
		double c_value[3];
		double dx;
	public:
		node1(void){};
		node1(double,int,double);
		node1(double,int,double,double,double,double);
		~node1(void);
		int GetLevel(void){return level;};
		double GetLoc(void){return loc;};
		double Getdx(void){return dx;};
		void GetVerts(double *);
		void GetCvals(double *);
		double GetI(void);
		void set_values(double (*)(double,void*),void*);
};
double node1::GetI(void){
	return (c_value[0]+c_value[1]+4*c_value[2])/6.;// Simpson
}
void node1::GetCvals(double *vals){
	for(int i=0;i<3;i++) vals[i]=c_value[i];
}
void node1::GetVerts(double *verts){
	verts[0]=loc; verts[1]=loc+dx; verts[2]=loc+.5*dx;
}
void node1::set_values(double (*fn)(double,void*),void* params){
	double points[3];
	GetVerts(points);
	for(int i=0;i<3;i++) c_value[i]=(*fn)(points[i],params);
}
node1::node1(double L,int lev,double DX,double value0,double value1,double value2){
	loc=L; level=lev; dx=DX/pow(2,level);
	c_value[0]=value0; c_value[1]=value1; c_value[2]=value2;
}
node1::node1(double lx,int lev,double DX){
	loc=lx; level=lev; dx=DX/pow(2,level);
}
node1::~node1(void){}

class location2{
	public:
		double x,y;
};

class intloc2{
	public:
		int x,y;
};

intloc2 fix2[4]={{0,0},{0,1},{1,1},{1,0}};

class node2{
	private:
		location2 loc;
		int level;
		double c_value[5];
		location2 diff;
	public:
		node2(void){};
		node2(double,double,int,location2);
		node2(location2,int,location2,double,double,double,double,double);
		~node2(void);
		int GetLevel(void){return level;};
		location2 GetLoc(void){return loc;};
		location2 Getdiff(void){return diff;};
		void GetVerts(location2 *);
		void GetCvals(double *);
		double GetSum(void);
		double GetVariance(double*);
		double GetI(void);
		void set_values(double (*)(double,double));
};

double node2::GetSum(void){
	double sum=c_value[0];
	for(int i=1;i<4;i++) sum+=c_value[i];
	return sum;
}
double node2::GetVariance(double *sum1){
	double sum2=pow(c_value[0],2); *sum1=c_value[0];
	for(int i=1;i<4;i++){
		*sum1+=c_value[i]; sum2+=pow(c_value[i],2);
	}
	return .25*(sum2-.25*pow(*sum1,2));
}
double node2::GetI(void){
	double sum=c_value[0];
	for(int i=1;i<4;i++) sum+=c_value[i];
	return (2*c_value[4]+sum)/6.;
}
void node2::GetCvals(double *vals){
	for(int i=0;i<5;i++) vals[i]=c_value[i];
}
void node2::GetVerts(location2 *verts){
	for(int i=0;i<4;i++){
		verts[i].x=loc.x+fix2[i].x*diff.x;
		verts[i].y=loc.y+fix2[i].y*diff.y;
	}
	verts[4].x=loc.x+.5*diff.x; verts[4].y=loc.y+.5*diff.y;
}
void node2::set_values(double (*fn)(double,double)){
	location2 points[5];
	GetVerts(points);
	for(int i=0;i<5;i++) c_value[i]=(*fn)(points[i].x,points[i].y);
}
node2::node2(location2 L,int lev,location2 DIFF,double value0,double value1,double value2,
	     double value3,double value4){
	loc=L; level=lev; double k=pow(2,level);
	diff.x=DIFF.x/k; diff.y=DIFF.y/k;
	c_value[0]=value0; c_value[1]=value1; c_value[2]=value2;
	c_value[3]=value3; c_value[4]=value4;
}
node2::node2(double lx,double ly,int lev,location2 DIFF){
	loc.x=lx; loc.y=ly; level=lev; double k=pow(2,level);
	diff.x=DIFF.x/k; diff.y=DIFF.y/k;
}
node2::~node2(void){}

class location3{
	public:
		double x,y,z;
};

class intloc3{
	public:
		int x,y,z;
};

intloc3 fix3[8]={{0,0,0},{0,1,0},{1,1,0},{1,0,0},{0,0,1},{0,1,1},{1,1,1},{1,0,1}};

class node3{
	private:
		location3 loc;
		int level;
		double c_value[9];
		location3 diff;
	public:
		node3(void){};
		node3(double,double,double,int,location3);
		node3(location3,int,location3,double,double,double,double,double,double,double,double,double);
		~node3(void);
		int GetLevel(void){return level;};
		location3 GetLoc(void){return loc;};
		location3 Getdiff(void){return diff;};
		void GetVerts(location3 *);
		void GetCvals(double *);
		double GetSum(void);
		double GetVariance(double*);
		double GetI(void);
		double GetIyb(void);
		double GetIxsq(void);
		double GetIysq(void);
		double GetIzsq(void);
		void set_values(double (*)(double,double,double));
};

double node3::GetSum(void){
	double sum=c_value[0];
	for(int i=1;i<8;i++) sum+=c_value[i];
	return sum;
}
double node3::GetVariance(double *sum1){
	double sum2=pow(c_value[0],2); *sum1=c_value[0];
	for(int i=1;i<8;i++){
		*sum1+=c_value[i]; sum2+=pow(c_value[i],2);
	}
	return .125*(sum2-.125*pow(*sum1,2));
}
double node3::GetI(void){// evaluation of rho
	double sum=c_value[0];
	for(int i=1;i<8;i++) sum+=c_value[i];
	return .125*(2*c_value[8]+.75*sum);
}
double node3::GetIyb(void){//for evaluation of rho*ybar
	double sum=0;
	for(int i=0;i<8;i++) sum+=c_value[i]*(loc.y+fix3[i].y*diff.y);
	return .125*(2*c_value[8]*(loc.y+0.5*diff.y)+.75*sum);
}
double node3::GetIxsq(void){//for evaluation of rho*xsqbar
	double sum=0;
	for(int i=0;i<8;i++) sum+=c_value[i]*pow(loc.x+fix3[i].x*diff.x,2);
	return .125*(2*c_value[8]*pow(loc.x+0.5*diff.x,2)+.75*sum);
}
double node3::GetIysq(void){//for evaluation of rho*ysqbar
	double sum=0;
	for(int i=0;i<8;i++) sum+=c_value[i]*pow(loc.y+fix3[i].y*diff.y,2);
	return .125*(2*c_value[8]*pow(loc.y+0.5*diff.y,2)+.75*sum);
}
double node3::GetIzsq(void){//for evaluation of rho*zsqbar
	double sum=0;
	for(int i=0;i<8;i++) sum+=c_value[i]*pow(loc.z+fix3[i].z*diff.z,2);
	return .125*(2*c_value[8]*pow(loc.z+0.5*diff.z,2)+.75*sum);
}
void node3::GetCvals(double *vals){
	for(int i=0;i<9;i++) vals[i]=c_value[i];
}
void node3::GetVerts(location3 *verts){
	for(int i=0;i<8;i++){
		verts[i].x=loc.x+fix3[i].x*diff.x;
		verts[i].y=loc.y+fix3[i].y*diff.y;
		verts[i].z=loc.z+fix3[i].z*diff.z;
	}
	verts[8].x=loc.x+.5*diff.x; verts[8].y=loc.y+.5*diff.y; verts[8].z=loc.z+.5*diff.z;
}
void node3::set_values(double (*fn)(double,double,double)){
	location3 points[9];
	GetVerts(points);
	for(int i=0;i<9;i++) c_value[i]=(*fn)(points[i].x,points[i].y,points[i].z);
}
node3::node3(location3 L,int lev,location3 DIFF,double value0,double value1,double value2,double value3,
	     double value4,double value5,double value6,double value7,double value8){
	loc=L; level=lev;
	double k=pow(2,level);
	diff.x=DIFF.x/k; diff.y=DIFF.y/k; diff.z=DIFF.z/k;
	c_value[0]=value0; c_value[1]=value1; c_value[2]=value2; c_value[3]=value3;
	c_value[4]=value4; c_value[5]=value5; c_value[6]=value6; c_value[7]=value7;
	c_value[8]=value8;
}
node3::node3(double lx,double ly,double lz,int lev,location3 DIFF){
	loc.x=lx; loc.y=ly; loc.z=lz; level=lev;
	double k=pow(2,level);
	diff.x=DIFF.x/k; diff.y=DIFF.y/k; diff.z=DIFF.z/k;
}
node3::~node3(void){}
