#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include <ctime>
#include <locale>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <thread>
#include <mutex>
#include <map>
#include <getopt.h>
#include <signal.h>

using namespace std;


vector<std::string> facs_names = { "Base R0 for virus", "Unreported", "Mobility multiply", "Daily Imported cases", "Infectious day mean", "Infectious day variance", "Cases reporting delay mean", "Case reporting delay variance", "Social Distancing effect",  "SD Intro day",  "SD Intro day var", "Start day" };
vector<double> facs_min = {1.5, 0.5, 0.5, 5,   1,    4,  5, 1, 0.6, 50,  10, 0  };
vector<double> facs_nom = {2.5, 0.8, 1,   10,  5,    6,  19, 4, 0.7, 75, 10, 1  };
vector<double> facs_max=  {5.0, 0.9, 1.5, 15,  10,   8,  25, 6, 0.9, 150,10, 30  };

void printfacsnames(ostream &outfile=cout) 
{
    for (const auto& e : facs_names) {
        outfile << setw(max(e.length()+2,14ul)) << e;
    }
    outfile << endl;
}

void printfacs(vector<double> &facs, ostream &outfile=cout) 
{
    for (int f=0; f<facs.size(); f++) {
        outfile << setprecision(9) << setw(max(facs_names[f].length()+2,14ul)) << facs[f];
    }
}


#define DAY0 "2019-12-30 12:00"
time_t day0=1577876400;
int day0_dow=3;
#define MAXDAYS 500

bool running=true;

int rcases[MAXDAYS]={0};
double mobility[MAXDAYS]={0};
double lastmobility=0;
int lastmobilityday=0;
int lastcaseday=0;
int weekendday=0;
double weekendeffect=0;


int strtday(string str)
{
    struct tm tm={};
    const char *cstr=str.c_str();
    if (strptime(cstr, "%Y-%m-%d", &tm)) {
        int d=((mktime(&tm) + (12*60*60))- day0) / (60*60*24);
        if (d<=0 || d>=MAXDAYS) {
            cout << "Day "<<cstr<<" read as day "<<d<<" doesn't make sense\n";
            cout << tm.tm_year<<" "<<tm.tm_mon<<" "<<tm.tm_mday<<"\n";
            exit(-1);
        }
        return d;
    } else {
        return -1;
    }
}
string daytstr(int d)
{
    char result[100];
    time_t t=(d * 60*60*24) + day0;
    strftime(result, 99, "%d/%m/%Y", gmtime(&t));

    return string(result);
}
int dow(int day)
{
    return (day+day0_dow)%7;
}

    

vector<string> split (const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;
    char chars[] = ",\n";

    while (getline (ss, item, delim)) {
        for (unsigned int i = 0; i < strlen(chars); ++i)
        {
            // you need include <algorithm> to use general algorithms like std::remove()
            item.erase (std::remove(item.begin(), item.end(), chars[i]), item.end());
        }
        result.push_back (item);
    }

    return result;
}




#define SQRT2PI 2.506628274631

double pdf ( double x, double m, double s )
{

    assert (s>0);

    double z = (x-m)/s;

    return exp(-0.5*z*z)/(SQRT2PI*s);
}

double cdf ( double x, double m, double s )
{

//  # Abramowitz & Stegun, 26.2.17
//  # absolute error less than 7.5e-8 for all x

    assert (s>0);

    double z = (x-m)/s;

    double t = 1.0/(1.0 + 0.2316419*abs(z));
    double y = t*(0.319381530
                 + t*(-0.356563782
                      + t*(1.781477937
                           + t*(-1.821255978
                                + t*1.330274429 ))));
    if( z > 0 ) {
        return 1.0 - pdf( z, 0, 1 )*y;
    } else {
        return pdf( z, 0, 1 )*y;
    }
}

double gaussian (double x, double mean, double var)
{
//    # adjust for end effects, assume numbers go 1-upwards
    if (x>(mean+(2*var))) { return 0; }
    double e=cdf(1, mean, var) + (1- cdf(mean+(2*var),mean,var));
    return pdf(x,mean,var) * (e+1);
}

mutex mtx;
map<double,pair<double,vector<double>>> bestruns;
vector<double> bestfacs;
double bestscore;
double bestcases;

double run (bool verbose, vector<double> &facs, bool record=false, ostream &outfile=cout)
{

    double R0_pop=facs[0];
    double unreported=facs[1]; // number of unreported cases 90%
    // Spain has ~50M 5% have antibodies, 2.5M, reported number 250k
    double mob_fac=facs[2];
    double importperday=facs[3];


    double infectday_m=facs[4];  // syptoms appear ~ day 5, 1 day before is most infectious
    double infectday_v=facs[5];


    double casedelay_m=facs[6];

    double casedelay_v=facs[7];
    double SD_effect=facs[8];
    double SD_introday=facs[9];
    double SD_introday_var=facs[10];

    //    day = Time::Piece->strptime("20191227", "%Y%m%d"); // day 0

    int startday=(int)facs[11];

    double casereports[100]={0};
    double infected[100]={0};

    double cases=0;
    double score=0;

    double oldcases=0;

    int future=40;
    

    vector<double> cpd(record?MAXDAYS:0);
    if (record) {
        fill(cpd.begin(), cpd.end(), 0);
    }
    
    if (verbose) outfile << "\"day\",\"date\",\"Estimated R0\", \"New infections\",\"best simulated cases\",\"Actual cases\", \"Actual delta\", \"Error\", \"Daily simulated delta\", \"Other simulated cases\"\n";
    
    for (int d=startday; d<lastmobilityday+future && d<lastcaseday+future; d++) {
        

            
        double R0=R0_pop;
        double SDe=1.0;
        /* This should be fixed*/
//        if (d-startday+1>=SD_introday) { SDe=SD_effect; }
        double f=cdf(d-startday+1, SD_introday, SD_introday_var);
        SDe=(f*SD_effect) + (1.0-f);
        
        if (d >= lastmobilityday) {
            R0=R0_pop * (((100.0+lastmobility) * mob_fac * SDe)/100.0);//R0;//_current;
        } else {
            if (mobility[d]) {
                R0=R0_pop * (((100.0+mobility[d]) * mob_fac * SDe)/100.0);
            } else {
                R0=R0_pop * (((100.0) * mob_fac * SDe)/100.0);
            }
        }


        int di=99-(d%100);
        infected[di]=0;
        casereports[di]=0;

        infected[di]+=importperday;

        double t=0;
        for( int i=1; i<100; i++) {
            double nc=R0 * infected[(di+i)%100] * gaussian(i, infectday_m, infectday_v);

            infected[di] += nc;

            casereports[di] += nc*(1-unreported);

            double cr;
            if (i<casedelay_m) {
                cr=casereports[(di+i)%100]*gaussian(i,casedelay_m, casedelay_v);
            } else {
                cr=casereports[(di+i)%100]*gaussian(casedelay_m,casedelay_m, casedelay_v);
            }
            if (dow(d)==weekendday) {
                cr*=weekendeffect;
            }
            cases+=cr;
            casereports[(di+i)%100]-=cr;
        }
        double err=0;
        if (d<lastcaseday && rcases[d]) {
            err=(rcases[d] - cases);
            if (d>lastcaseday-7) {
                err*=2;
            } else {
                err=min(abs(err),cases/2.0);
            }
            
            score+=(err*err)/((lastcaseday+1)-d);
        }
        if (record) cpd[d]=cases;
        if (verbose) {
            outfile << d-startday+1<<", "<<daytstr(d)<<", "<<R0<<", "<<infected[di]<<", "<<cases<<", "
                    <<((d<lastcaseday) && rcases[d]?to_string(rcases[d]):" ")<<", "
                 <<(rcases[d]&&rcases[d-1]?to_string(rcases[d]-rcases[d-1]):" ")
                 <<", "<<err<<", "<<(cases-oldcases);
            for (auto itr = bestruns.crbegin(); itr != bestruns.crend(); ++itr) { 
                if (itr->second.second.size()>d) {
                    if (itr->second.second[d]>(cases*1.5)) outfile <<", ";
                    else outfile << ", " << itr->second.second[d];
                }
                
            } 
            outfile <<"\n";
        }
        
        oldcases=cases;
        // day +1
    }

    if (record && (abs(bestscore-score)/bestscore < 1)) {
        if (score < bestscore) {
            bestcases=cases;
//            bestruns.clear();
        } else {   
            unique_lock<mutex> lck(mtx);
            double s= -(abs(bestcases-cases)/abs(bestscore-score));
            bestruns.insert(make_pair(s,make_pair(score,cpd)));
            while (bestruns.size()>10) {
                bestruns.erase(--bestruns.end());
            }
        }
    }
  
    
    
    if (verbose) {
        outfile<< "Score: "<< score<<"   , , , , , , , , ";
        for (auto itr = bestruns.crbegin(); itr != bestruns.crend(); ++itr) { 
            outfile << ", \"score: " << itr->second.first << " \"";
        } 
        outfile <<"\n";
    }
    
    return score;
}

double findlocalmin (vector<double> &facs, double best)
{

    int dinit=1000;
    for (int i=0;i<10;i++) {
        int changes=0;
        for (int f=0; f<facs.size(); f++) {
            if (facs_min[f]!=facs_max[f]) {
                int d=dinit;
                int found=0;
                while (abs(d)<=10000 && found<10) {
                    double old=facs[f];
                    facs[f] += facs[f]*(1.0/(double)d);
                    if (facs[f]<facs_min[f]) facs[f]=facs_min[f];
                    if (facs[f]>facs_max[f]) facs[f]=facs_max[f];
                    
                    double s=run(0, facs);
                    if (s < best) {
//                        print "round i D=d Change=".(best-s)." ".facs_names[f]."\n";
                        changes+=(best-s);
                        best=s;
                        found+=1;
                        if (found>=5) {
                            d/=2;
                            found=0;
                            if (abs(d)<=2) {
                                return best;
                            }
                        }
                    } else {
                        facs[f]=old;
                        if ((d>0) && found==0) {
                            d*=-1;
                        } else {
                            found=0;
                            d=abs(d*10);
                            if (facs[f]==facs_max[f]) d*=-1;
                        }
                    }
                }
            }
        }
//        print "round i Changes=changes\n";
        if (changes<500) break;
//        print "best\n";
        dinit*=5;
    }
    return best;
}


double better(bool verbose, vector<double> &test)
{
    double bests=run(0, test);
    double old;
    if (verbose) printfacsnames();
    do {
        old=bests;
        bests =findlocalmin(test, bests);
        if (verbose) {
            printfacs(test);
            cout << " Score= "<<bests<<" ("<<(bests-old)<<")\n";
        }
    } while (old!=bests);

    if (verbose) {
        bests =run(true, test);
//        print join(' ',@facs_names),"\n";
//        print join(' ',@facs),"\n";
        printfacsnames();
        printfacs(test);
        cout << "\n Score= "<<bests<<"\n";
    }
    return bests;
}

int testid=0;
void runandtests(bool verbose=true, int runs=10000, int ignoreAbove=100)
{
    double temp=1;
    
    vector<double> test=facs_nom;
//    @test=@facs;
    while (testid < runs && running) {// && temp>0.1) {
        {
            unique_lock<mutex> lck(mtx);
            testid+=1;
            if (verbose) {
                printf("\r%7d ",testid);
                cout << flush;
            }
            temp=(1.0 - (double)testid/((double)runs*1.1) );
            test=facs_nom;
            vector<double> current=bestfacs;
            for (int i=0; i<facs_nom.size(); i++) {
//                if (rand(100000) > testid/2) {
                if (facs_min[i]==facs_max[i]) {
                    test[i]=facs_max[i];
                } else {
                    double max=current[i]+((facs_max[i]-current[i]) * temp);
                    double min=current[i]-((current[i]-facs_min[i]) * temp);
//                    max=facs_max[i];
//                    min=facs_min[i];
//                    test[i] = rand(max-min)+min;
                    test[i] = min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max-min)));
                }
//                }
            }
//            printfacs(test);

        }
        double locals=run(0, test);
        {
            if (locals < bestscore*ignoreAbove) {
//                double s=findlocalmin(test, locals);
                double s=better(false, test);
                run(false, test, true); // remember the scores !
                {
                    unique_lock<mutex> lck(mtx);
                    
                    if (verbose) {
                        printfacs(test);
                         cout << " Score="<<setw(12)<<s<<" ("<<(locals/s)<<") ";
                    }
                    

                    if (s<bestscore) {

                        if (verbose) cout<< "(*)";
                        bestscore=s;
                        bestfacs=test;

//                        temp*=0.95;
                    }
                }
                if (verbose) cout<< "\n";
            }
        }
    }
}

void intHandler(int dummy) {
    running=false;
}

vector<double> argvtv(int start, int argc, char* argv[]) 
{
    vector<double> result;
    for (;start<argc;start++) {
        result.push_back(stod(argv[start]));
    }
    return result;
}

int main(int argc, char *argv[])
{
    std::srand(std::time(0));

    signal(SIGINT, intHandler);

    cout << setprecision(9);
    
    {
        struct tm tm={};
        const char *d=DAY0;
        strptime(d, "%Y-%m-%d %H:%M", &tm);
        day0=mktime(&tm);
        
        day0_dow=tm.tm_wday;
    }
    

    int ignoreAbove=100;
    int stopday=MAXDAYS; // stop reading cases after this day
    bool cont=false;
    bool facsprovided=false;
    ofstream outfile;
    int runs=10000;
    bool verbose=true;
    string country="unknown";
    int c;
    do {
        while ((c=getopt(argc,argv,"co:n:qs:i:")) != -1)
        {
            switch (c) {
                case 'c':
                    cont=true;
                    break;
                case 'o':
                    outfile.open(optarg);
                    outfile << setprecision(9);
                    break;
                case 'n':
                    runs=stoi(optarg);
                    break;
                case 'q':
                    verbose=false;
                    break;
                case 's':
//                    stopday=strtday(optarg);
                    stopday=stoi(optarg);
                    break;
                case 'i':
                    ignoreAbove=stoi(optarg);
                    break;
            }
        }
        if (country=="unknown") {
            if (optind < argc) {
                country=argv[optind++];
                continue;
            } else {
                cout << "Usage :" << argv[0] << " country [-c][-o file.csv][-n runs][-q] [factors]\n"; exit(-1);
            }
        } else {
            break;
        }
    } while (true);
    
    if (optind+1 < argc) {
        facsprovided=true;
        facs_nom=argvtv(optind,argc,argv);
    }
    
/*
    int casesperday[8]={0};
    ifstream myfile (country+".csv");
    int oldcases=0;
    if (myfile.is_open()) {
        string line;
        while(getline(myfile,line))
        {
            vector<string> fields = split(line,' ');
            int d=strtday(fields[0]);
            
            if (d>=0) {
                int cases=stoi(fields[1]);
                rcases[d]=cases;
                casesperday[dow(d)]+=cases-oldcases;
                oldcases=cases;
            }
        }
        myfile.close();
    } else {
        cout << "Unable to open "+country+".csv\n"; exit(-1);
    }
*/
    int casesperday[8]={0};
    ifstream myfile ("owid-covid-data.csv");
    int oldcases=0;
    int found=0;
    if (myfile.is_open()) {
        string line;
        while(getline(myfile,line))
        {
            vector<string> fields = split(line,',');
            if (fields[2]==country) {
                int d=strtday(fields[3]);
                if (d>=0 && d<MAXDAYS) {
                    int cases=stoi(fields[4]);
                    rcases[d]=cases;
                    casesperday[dow(d)]+=cases-oldcases;
                    oldcases=cases;
                    lastcaseday=d;
                    found++;
                }
                
            }
            
        }
        myfile.close();
    } else {
        cout << "Unable to open owid-covid-data.csv\n"; exit(-1);
    }
    if (!found) {
        cout << "Unable to find country "<<country<<" in owid-covid-data.csv\n"; exit(-1);
    }

    if (lastcaseday>stopday) lastcaseday=stopday;
    
    
    int total=0;
    for (int d=0;d<7;d++) {
        if (casesperday[d] < casesperday[weekendday]) weekendday=d;
        total+=casesperday[d];
    }
    weekendeffect=(double)(casesperday[weekendday])/((double)total/7.0);
    
/*    ifstream data ("mobility-"+country+".csv");
    if (data.is_open()) {
        string line;
        while (getline(data,line)) {
            vector<string> fields = split(line,',');
            
            int d=strtday(fields[4]);
            if (d>=0) {
                double m=(double)(stoi(fields[5])+stoi(fields[6])+stoi(fields[7])+stoi(fields[8])+stoi(fields[9])-stoi(fields[10]))/6.0;
                mobility[d]=m;
                lastmobility=((lastmobility*2.0)+m)/3.0;
                lastmobilityday=d;
            }
        }
        data.close();
    } else {
        cout<< "Cant find mobility-"+country+".csv\n"; exit(-1);
    }
*/
    found=0;
    ifstream data ("Global_Mobility_Report.csv");
    if (data.is_open()) {
        string line;
        while (getline(data,line)) {
            vector<string> fields = split(line,',');
            if (fields[1]==country && fields[2]=="") {
                int d=strtday(fields[7]);
                if (d>=0) {
                    double m=(double)(stoi(fields[8])+stoi(fields[9])+stoi(fields[10])+stoi(fields[11])+stoi(fields[12])-stoi(fields[13]))/6.0;
                    mobility[d]=m;
                    lastmobility=((lastmobility*2.0)+m)/3.0;
                    lastmobilityday=d;
                    found++;
                }
            }
        }
        data.close();
    } else {
        cout<< "Cant find Global_Mobility_Report.csv\n"; exit(-1);
    }
    if (!found) {
        cout << "Unable to find country "<<country<<" in Global_Mobility_Report.csv\n"; exit(-1);
    }

    vector<double> cpd(MAXDAYS);
    fill(cpd.begin(), cpd.end(), 0);
    for (int i=0; i<10; i++) 
    {
        bestruns.insert(make_pair(0,make_pair(0,cpd)));
    }
    

    if (!cont && facsprovided) {
        if (outfile.is_open()) {
            run(true,facs_nom,false,outfile);
            outfile.close();
        } else {
            run(true,facs_nom);
        }
        exit(-1);
    }
    bestfacs  =facs_nom;
    bestscore =run(false, facs_nom);


    if (verbose) {
        cout << "        "; // allow for test number
        printfacsnames();
    }
    

    testid=0;
        
    vector<thread> threads;
    for (int t=0;t<std::thread::hardware_concurrency();t++) {
        threads.push_back(thread(runandtests,verbose,runs, ignoreAbove));
    }
    while (!threads.empty()) {
        threads.back().join();
        threads.pop_back();
    }

    double score=better(verbose, bestfacs);

    if (outfile.is_open()) {
        bestscore =run(true, bestfacs, false, outfile);
        printfacsnames(outfile);
        printfacs(bestfacs,outfile);
        outfile.close();
    }
    
}


