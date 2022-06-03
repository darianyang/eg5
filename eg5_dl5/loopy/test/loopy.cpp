#include"../fasp/source.h"
Rcs RCS("model.dir");
Tres TRES;

main(int argc,char *argv[])
{
Strhandler cc;
Segen segen; 
Chiangle chiangle;
Pdb pdb;
char name[100],loopname[1000];
char *seqn=new char[1000];
name[0]='\0';
loopname[0]='\0';
seqn[0]='\0';
int kk=0;
TRES.flag=102;
TRES.smt=1;
if(argc==1) {
segen.printhelp();exit(0);

}
int n;
for(n=1;n<argc;n++)
{
  if(strstr(argv[n],"-o="))  //loops to be predicted
  {
    segen.start=atoi(argv[n]+3);
    char *poi=strstr(argv[n]+3,"-")+1;
    if(poi==0) {
       cerr<<"error in command specifying start and end residues"<<endl;
       segen.printhelp();
       exit(0);
    }
    segen.end=atoi(poi);
  }
  else if(strstr(argv[n],"-h")) {
    segen.printhelp();
  }
  else if(strstr(argv[n],"-i="))  //number of initials 
  {
    segen.arbt=atoi(argv[n]+3);
    if(segen.arbt==0) {
	cerr<<"number of initials should not be zero:"<<argv[n]<<endl;
        exit(0);
    }
    segen.part=segen.arbt/10;
    segen.part0=segen.arbt/2;
    if(segen.part<10 ) segen.part = segen.part0;
  }
  else if(strstr(argv[n],"-t=")) //topology
  {
    
    if(argv[n][3]=='a') TRES.flag=100;
    else if(argv[n][3]=='b')  TRES.flag=102;    
    else if(argv[n][3]=='A')  TRES.flag=200;
    else if(argv[n][3]=='B')  TRES.flag=202; 
    
    //TRES.flag=atoi(argv[n]+3);
  }
  else if(strstr(argv[n],"-n=")) //number of outputs
  {
    segen.revs=atof(argv[n]+3);
    if(segen.revs==0) {
	cerr<<"number of output should not be zero:"<<argv[n]<<endl;
        exit(0);
    }
  }
  else if(strstr(argv[n],"-s=")) { //secondary structure
    segen.secd=argv[n][3];
    if(segen.secd=='c') segen.secd='-';
    if(segen.secd!='-'&&segen.secd!='h'&&segen.secd!='e') {
	cerr<<"secondary structure could only be e or h or -"<<endl;
        exit(0);
    }
    
  }
  else if(strstr(argv[n],"-c=")) { //chain id
    segen.cid=argv[n][3];
  }
  else if(strstr(argv[n],"-r=")) {
    	strcpy(seqn,argv[n]+3);    //sequence
	seqn=cc.clearendchar(seqn," \n");
	if(seqn==0) {
		cerr<<argv[n]<<" is not a sequence"<<endl;	
		exit(0);	
	} 
        seqn=cc.uppercase(seqn);
  }
  else if(strstr(argv[n],"-f")) {
	segen.fapr=1;
  }
  else if(strstr(argv[n],"-")) {
    continue;
  }
  else
  {
    kk++;
    if(kk==1) strcpy(name,argv[n]);         //protein name
    if(kk==2) {
		strcpy(loopname,argv[n]);
		segen.databaseonly=1;
		
    }
  }
}


TRES.read("tres");
TRES.readmore("back_small_rotamer","side_small_rotamer");

if(strlen(loopname)!=0) {
	chiangle.read(loopname);
	segen.chiangle=&chiangle;
}

seqn=cc.clearendchar(seqn,"\"");

if(seqn) {
if(strlen(seqn)!=segen.end-segen.start+1) {
	cerr<<"number of residues in loop from "<<segen.start<<" to "\
	<<segen.end<<" does not match the sequnece given"<<endl;
	exit(0);
}
for(int ii=0;ii<strlen(seqn);ii++) {
	if(!TRES[seqn[ii]]) {
		cerr<<seqn[ii]<<" is not a standard protein residue"<<endl;
		exit(0);
	}
}
}


Tres *tt;
Tatm *aa;
for(tt=&TRES;tt;tt=tt->next)
for(aa=tt->tatm;aa;aa=aa->next)
{
//aa->eng->radius=2.0;
cerr<<tt->name<<" "<<aa->name<<" "<<aa->eng->radius<<"  "<<aa->eng->charge<<endl;
}

pdb.read(name,'0');
if(pdb.chn==0) exit(0);
if(pdb.number>1&&segen.cid=='1') {
  cerr<<"please specify the chain id on which loop is predicted.."<<endl;
  exit(0);
}
segen.pdb=&pdb;

if(segen.cid=='1') segen.cid=pdb.chn->id;
Chn *chn=pdb[segen.cid];
if(chn==0) {
   cerr<<"chain: "<<segen.cid<<" does not exist!"<<endl;
   exit(0); 
}


//add missing residues
if(seqn&&seqn[0]!='\0') chn->addresiduesonly(segen.start,segen.end,seqn);

chn->write("re");
Res *r;

r=chn->isres(segen.start-chn->start);
if(r==0) {
  cerr<<"residue: "<<segen.start<<" in chain: "<<segen.cid<<" does not exist!"<<endl;
  exit(0);
}
segen.start=r->id0;

r=chn->isres(segen.end-chn->start);
if(r==0) {
  cerr<<"residue: "<<segen.end<<" in chain: "<<segen.cid<<" does not exist!"<<endl;
  exit(0);
}
segen.end=r->id0;
 
segen.ready();
float **xyz_all=segen.predt0();

n=0;
while(xyz_all&&xyz_all[n])n++;

if(n>=1)
{
  int ii=0;
  FILE *fp=stdout;
  for(int i=0;i<n;i++)
  {
    if(ii>=segen.revs) break;
    segen.segment->chn->transfer(xyz_all[i],1);
    float d=segen.rmsd(2);
    if(d<0.001) continue;
    char ch1=segen.cid;
    if(ch1==' ') ch1='0';
    fprintf(fp,"vw: %2i %s %c %5i %5i %8.3f\n",i,name,ch1,segen.start,segen.end,d);
    ii++;
  }
  fclose(fp);
}

xyz_all=cc.floatdel(xyz_all);
if(seqn) delete [] seqn;
exit(0);
}
