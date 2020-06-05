#include <iostream>
#include <sstream>
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

using namespace std;


vector<std::string> facs_names = { "Base R0 for virus", "Unreported", "Mobility multiply", "Daily Imported cases", "Infectious day mean", "Infectious day variance", "Cases reporting delay mean", "Case reporting delay variance", "Social Distancing effect",  "SD Intro day", "Start day" };
vector<double> facs_min = {1.5, 0.1, 0.5, 1,   5,   5, 15, 3, 0.5, 50,  0  };
vector<double> facs_nom = {2.5, 0.9, 1,   10,  5,   5, 19, 4, 0.8, 100, 1  };
vector<double> facs_max=  {5.0, 1.0, 2.5, 20,  5,   5, 25, 6, 1.2, 150, 30  };

void printfacsnames() 
{
    for (const auto& e : facs_names) {
        cout << setw(max(e.length()+2,14ul)) << e;
    }
    cout << endl;
}

void printfacs(vector<double> &facs) 
{
    for (int f=0; f<facs.size(); f++) {
        cout << setprecision(9) << setw(max(facs_names[f].length()+2,14ul)) << facs[f];
    }
}


#define DAY0 "2020-01-01 12:00"
time_t day0=1577876400;
int day0_dow=3;
#define MAXDAYS 500


int rcases[MAXDAYS]={0};
double mobility[MAXDAYS]={0};
double lastmobility=0;
int lastmobilityday=0;
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
    char chars[] = ", \t\n";

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
map<double,vector<double>> bestruns;

double run (bool verbose, vector<double> &facs, bool record=false)
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

    //    day = Time::Piece->strptime("20191227", "%Y%m%d"); // day 0

    int startday=(int)facs[10];

    double casereports[100]={0};
    double infected[100]={0};

    double cases=0;
    double score=0;

    double oldcases=0;

    vector<double> cpd(record?MAXDAYS:0);
    if (record) {
        fill(cpd.begin(), cpd.end(), 0);
    }
    
    if (verbose) cout << "\"day\",\"date\",\"Estimated R0\", \"New infections\",\"best simulated cases\",\"Actual cases\", \"Actual delta\", \"Error\", \"Daily simulated delta\", \"Other simulated cases\"\n";
    
    for (int d=startday; d<lastmobilityday+25; d++) {
        

            
        double R0=R0_pop;
        double SDe=1.0;
        /* This should be fixed*/
        if (d-startday+1>=SD_introday) { SDe=SD_effect; }
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
        if (rcases[d]) {
            err=(rcases[d] - cases);
            score+=(err*err)/(200.0-(d-startday+1));
        }
        if (record) cpd[d]=cases;
        if (verbose) {
            cout << d-startday+1<<", "<<daytstr(d)<<", "<<R0<<", "<<infected[di]<<", "<<cases<<", "
                 <<(rcases[d]?to_string(rcases[d]):" ")<<", "
                 <<(rcases[d]&&rcases[d-1]?to_string(rcases[d]-rcases[d-1]):" ")
                 <<", "<<err<<", "<<(cases-oldcases);
            for (auto itr = bestruns.crbegin(); itr != bestruns.crend(); ++itr) { 
                if (itr->second.size()>d) cout << ", " << itr->second[d]; 
            } 
            cout <<"\n";
        }
        
        oldcases=cases;
        // day +1
    }

    if (record) {
        unique_lock<mutex> lck(mtx);
        bestruns.insert(make_pair(score,cpd));
        while (bestruns.size()>10) {
            bestruns.erase(--bestruns.end());
        }
    }
  
    
    
    if (verbose) {
        cout<< "Score: "<< score<<"   , , , , , , , , ";
        for (auto itr = bestruns.crbegin(); itr != bestruns.crend(); ++itr) { 
            cout << ", \"score: " << itr->first << " \"";
        } 
        cout <<"\n";
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
    double bestscore=run(0, test);
    double old;
    if (verbose) printfacsnames();
    do {
        old=bestscore;
        bestscore =findlocalmin(test, bestscore);
        if (verbose) {
            printfacs(test);
            cout << " Score= "<<bestscore<<" ("<<(bestscore-old)<<")\n";
        }
    } while (old!=bestscore);

    if (verbose) {
        bestscore =run(true, test);
//        print join(' ',@facs_names),"\n";
//        print join(' ',@facs),"\n";
        printfacsnames();
        printfacs(test);
        cout << "\n Score= "<<bestscore<<"\n";
    }
    return bestscore;
}

vector<double> bestfacs;
double bestscore;
int testid=0;
void runandtests(bool verbose=true)
{
    double temp=1;
    
    vector<double> test=facs_nom;
//    @test=@facs;
    while (testid < 95000) {// && temp>0.1) {
        {
            unique_lock<mutex> lck(mtx);
            testid+=1;
            if (verbose) {
                printf("\r%7d ",testid);
                cout << flush;
            }
            temp=(1 - testid/100000 );
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
            if (locals < bestscore*100) {
//                double s=findlocalmin(test, locals);
                double s=better(false, test);
                run(false, test, true); // remember the scores !
                {
                    unique_lock<mutex> lck(mtx);
                    
                    if (verbose) {
                        printfacs(test);
                        cout << " Score="<<s<<" ";
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

    cout << setprecision(9);
    
    {
        struct tm tm={};
        const char *d=DAY0;
        strptime(d, "%Y-%m-%d %H:%M", &tm);
        day0=mktime(&tm);
        
        day0_dow=tm.tm_wday;
    }
    
    if (argc ==1) {
        cout << "Usage :" << argv[0] << " country [-c][-b] [factors]\n"; exit(-1);
    }
    
    string country=argv[1];

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

    int total=0;
    for (int d=0;d<7;d++) {
        if (casesperday[d] < casesperday[weekendday]) weekendday=d;
        total+=casesperday[d];
    }
    weekendeffect=(double)(casesperday[weekendday])/((double)total/7.0);
    
    ifstream data ("mobility-"+country+".csv");
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

    bool verbose=true;

    if (argc > 3 && strcmp(argv[2], "-b")==0) {
        
        vector<double> test=argvtv(3,argc,argv);
        better(true, test);
        exit(-1);
    }

    if (argc > 3 && strcmp(argv[2], "-c")!=0) {
        vector<double> test=argvtv(2,argc,argv);
        run(true, test);
        exit(-1);
    }

    
    if (argc > 3 && strcmp(argv[2], "-c")==0) {
        facs_nom=argvtv(3,argc,argv);
        verbose=false;
//        cout<< "Continuing...\n";
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
        threads.push_back(thread(runandtests,verbose));
    }
    while (!threads.empty()) {
        threads.back().join();
        threads.pop_back();
    }

    double score=better(verbose, bestfacs);

    if (!verbose) {
        bestscore =run(true, bestfacs);
        printfacsnames();
        printfacs(bestfacs);
    }
}




//1.88655952464931 1.29484178509565 10.711529000126 5.72632099401432 5.37355998352172 12.9839263340605 0.903301995852417 0.700699272597673 17.9150521022227  Score= 1008967.16566595
//./covid.pl france 1.62315429220329 1.22322391786631 14.6589441000495 3.35628100970738 3.81210295741958 13.7014759097719 3.34789042931359 1.02949830353047 23.6918429689992 score = 983152.076326347

//./covid.pl france 1.62315429220329 0.9 1.22322391786631 14.6589441000495 3.35628100970738 3.81210295741958 13.7014759097719 3.34789042931359 1.02949830353047 75.0 23.6918429689992 score = 983152.076326347
//./covid.pl france 1.62315429220329 0.9 1.22322391786631 14.6560124578189 3.35628100970738 3.81210295741958 13.7055866265469 3.34153194906842 1.02918948492333 75.0 23.6918429689992  Score= 985718.6959901




