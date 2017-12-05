#include<stdio.h>
#include<stdlib.h>
#include<cstdlib>
#include<time.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<iomanip>
#include<vector>
#include<math.h>
#include<algorithm>

using namespace std;

struct individual {
	int num;
	int x;
	int y;
	int demand;
};

vector<individual> individuals; //requirements
vector< vector<int> > population;
int individual_num = 250;
int population_size = 200;
int capacity = 500;

double best_distance = 99999;
double crossoverprob = 0.7;
double mutationprob = 0.02;
double insertionprob = 0.02;
double greedyprob = 0.1;
vector<int> cur_chromosome;

vector<int> best_chromosome;
double mintrucksolution = 99999;
//vector<int> minroute; //minroute for current truck
double delta = 0;
double T = 2; //initial temperature
double a = 0.95;
double k = 0.95;

void readfile(string filename)
{
	string lineA;
	string coorsec = "NODE_COORD_SECTION";
	string demasec = "DEMAND_SECTION";
	ifstream fileIn;
	int flag = 0;
	int row1 = 0, row2 = 0;


	//open file and error check
	fileIn.open(filename.c_str());
	if (fileIn.fail()) {
		cerr << " * the file you are trying to access cannot be found or opened";
		exit(1);
	}

	//initialize vector
	individual temp;
	for (int i = 0; i <= individual_num; i++) {
		individuals.push_back(temp);
	}


	//Read data file to two arrays
	while (fileIn.good()) {
		while (getline(fileIn, lineA)) {

			//cout << lineA << endl;

			// if found scanning coordinates section, set flag to 1 and read next line
			if (lineA.find(coorsec) < lineA.length()) { flag = 1; continue; }
			// if found scanning demands section, set flag to 2 and read next line
			if (lineA.find(demasec) < lineA.length()) { flag = 2; continue; }

			//when flag equals to 1 read data to coordinates array
			if (flag == 1) {
				//cout << "flag ==1" << endl;				
				istringstream stream(lineA);          //set separators to the line string 
				int col1 = 0;                         //so as to recognize each number

													  //put number and coordinates of individuls to the vector
				stream >> individuals[row1].num >> individuals[row1].x >> individuals[row1].y;
				//cout << individuals[row1].num << individuals[row1].x << individuals[row1].y;
				row1++;                               // go to next row of the array
			}

			//when flag equals to 2 read data to demands 
			if (flag == 2) {
				//cout << "flag ==2" << endl;		
				istringstream stream(lineA);
				stream >> individuals[row2].num >> individuals[row2].demand;
				//cout << individuals[row2].demand;				
				row2++;
			}
		}
	}
}

void initialize() {

	vector<int> counter;
	for (int j = 0; j < individual_num; j++)      //counter has to be 1 element larger than chromosome, because counter[0] is not in use
		counter.push_back(0);                     //and we only want generate random number between 1 to 249, so it must have counter[249]

												  //generate random number from 1 to 249 (client reference number stored in vector 'individuals'(location 0 is depot) and push them into chromosome
	for (int cnt = 1; cnt < individual_num; cnt++) {

		//int r = rand() % (MAX - MIN + 1) + MIN;
		int r = rand() % (individual_num - 1) + 1;

		while (++counter[r] > 1)
			r = rand() % (individual_num - 1) + 1;

		cur_chromosome.push_back(r);
	}

}

vector<int> insert_depot(vector<int> chromosome) {
	vector<int> route_map;
	route_map.push_back(0);
	int capacity_used = 0;

	for (int i = 0; i <chromosome.size(); i++) {

		if (capacity_used + individuals[chromosome[i]].demand < capacity) {
			route_map.push_back(chromosome[i]);
			capacity_used = capacity_used + individuals[chromosome[i]].demand;
		}
		else {
			route_map.push_back(0);
			route_map.push_back(chromosome[i]);
			capacity_used = individuals[chromosome[i]].demand;
		}
	}

	route_map.push_back(0);
	return route_map;
}

double calculate_distance(vector<int> chromosome) {
	vector<int> route_map = insert_depot(chromosome);
	double distance = 0;
	for (int i = 1; i < route_map.size(); i++) {
		double xdistance = individuals[route_map[i]].x - individuals[route_map[i - 1]].x;
		double ydistance = individuals[route_map[i]].y - individuals[route_map[i - 1]].y;
		distance = distance + sqrt((pow(xdistance, 2) + pow(ydistance, 2)));
	}
	return distance;
}

string num2str(int i)
{
         stringstream ss;
         ss<<i;
         return ss.str();
 }
int main() {

	clock_t start, finish;
	double totaltime;
	start = clock();




	int max_generation=100000;
	cout << "Please input maximum generation" << endl;
	//cin >> max_generation;

	// set current time as random seed
	srand((int)time(0) + rand());

	//read file and store in vectors
	readfile("fruitybun250.vrp");

	
	ofstream evolutioncurve;
	evolutioncurve.open("evolutioncurve.txt",ios::trunc);

	ofstream currentbest;
	currentbest.open("currentbest.txt",ios::trunc);
	currentbest<<"login fr16481 1347322"<<endl;
	currentbest<<"name Feifei Rong"<<endl;
	currentbest<<"algorithm SA with mutation, inverse and insertion"<<endl;
	
	ofstream plot_cvrp;
	plot_cvrp.open("plot-cvrp",ios::trunc);
	plot_cvrp << "set term png"<<endl;
	plot_cvrp <<"set output 'solution.png'"<<endl;
	plot_cvrp << "set size 1, 1"<<endl;
	plot_cvrp << "set grid"<<endl;
	plot_cvrp << "plot 'evolutioncurve.txt' using 1:2 with linespoints lt 3 lw 2 pt 7 ps .5" <<endl;
	plot_cvrp << "pause -1"<<endl;

	ofstream plot_route;
	plot_route.open("plot-route",ios::trunc);
	plot_route << "set term png"<<endl;
	plot_route << "set output 'route.png'"<<endl;

	plot_route << "set size 1,1"<<endl;
	plot_route <<"set pointsize 1.2"<<endl;
	plot_route <<"set xrange [-120:120]"<<endl;
	plot_route <<"set yrange [-120:120]"<<endl;
	plot_route <<"set style data linespoints"<<endl;
	plot_route <<" plot 'truck1.txt' using 2:3 t 'truck 1', 'truck2.txt' using 2:3 t 'Truck 2', 'truck3.txt' using 2:3 t 'Truck 3', 'truck4.txt' using 2:3 t 'truck 4', 'truck5.txt' using 2:3 t 'truck 5','truck6.txt' using 2:3 t 'truck 6','truck7.txt' using 2:3 t 'Truck 7','truck8.txt' using 2:3 t 'truck 8','truck9.txt' using 2:3 t 'Truck 9','truck10.txt' using 2:3 t 'truck 10','truck11.txt' using 2:3 t 'truck 11','truck12.txt' using 2:3 t 'truck 12','truck13.txt' using 2:3 t 'truck 13','truck14.txt' using 2:3 t 'truck 14','truck15.txt' using 2:3 t 'truck 15','truck16.txt' using 2:3 t 'Truck 16','truck17.txt' using 2:3 t 'Truck 17','truck18.txt' using 2:3 t 'truck 18','truck19.txt' using 2:3 t 'Truck 19','truck20.txt' using 2:3 t 'truck 20','truck21.txt' using 2:3 t 'truck 21','truck22.txt' using 2:3 t 'truck 22','truck23.txt' using 2:3 t 'truck 23','truck24.txt' using 2:3 t 'truck 24','truck25.txt' using 2:3 t 'Truck 25','truck26.txt' using 2:3 t 'truck 26' "<<endl;
	plot_route << "pause -1"<<endl;

	initialize();

	double cur_distance;
	int generation = 1;
	while (generation < max_generation) {
		
		//mutation 
		//if(rand() < 0.5){
		for(int j=0;j<250;j++){

			//old distance
			vector<int> old_chromosome = cur_chromosome;
	
			double distance = calculate_distance(cur_chromosome);
			double old_distance = distance;
			if (distance < best_distance)
				best_distance = distance;
			double r= rand()%10000*0.0001;
			vector<int> new_chromosome = cur_chromosome;

			if(r <= 0.25){
				//start mutation
				int index1 = rand() % (individual_num - 1);
				int index2 = rand() % (individual_num - 1);

				while (index1 == index2)
					index2 = rand() % (individual_num - 1);

				int temp = new_chromosome[index1];
				new_chromosome[index1] = new_chromosome[index2];
				new_chromosome[index2] = temp;
			}

			else if(r >= 0.5){
				int index1 = rand() % (individual_num - 2)+1;
				int index2 = rand() % (individual_num - 2)+1;

				while (index1 == index2)
					index2 = rand() % (individual_num-2)+1;

				while (index1 > index2) {
					int temp = index2;
					index2 = index1;
					index1 = temp;
				}

				int pickout = cur_chromosome[index2];
				vector<int> new_chromosome = cur_chromosome;
				new_chromosome.erase(new_chromosome.begin() + (index2));
				new_chromosome.insert(new_chromosome.begin() + (index1),pickout);
			}

			else{
				int index1 = rand() % (individual_num - 1);
				int index2 = rand() % (individual_num - 1);

				while (index1 == index2)
					index2 = rand() % (individual_num - 1);

				//index1 is start gene, index2 is end gene
				while (index1 > index2) {
					int temp = index1;
					index1 = index2;
					index2 = temp;
				}
				
		
				//give inverse part to inversetemp
				vector<int> inversetemp;
				for (int j = index1; j <= index2; j++) {
					inversetemp.push_back(new_chromosome[j]);
				}
				//inverse in newchromosome
				int k = inversetemp.size() - 1;
				for (int j = index1; j <= index2; j++) {
					new_chromosome[j] = inversetemp[k];
					k = k - 1;
				}
			}		
			
			distance = calculate_distance(new_chromosome);
			double new_distance = distance;

			if (distance < best_distance){
				best_distance = distance;
				best_chromosome.clear();
				best_chromosome = new_chromosome;
			}
			//apply SA
			delta = new_distance - old_distance;
			if (delta > 0) { //if has no progress
				double rnd = rand() % 10000 * 0.00001;
				double p = exp(-delta / (k*T));//possibility of accept new population is p.
				if (rnd < p) { //when random number lie within p, update chromosome
				cur_chromosome = new_chromosome;
				}
			}
			else
			cur_chromosome = new_chromosome;

		}
		cur_distance = calculate_distance(cur_chromosome);
		
		evolutioncurve << generation << " " << cur_distance << endl;
		T = T*a;

		generation++;
		cout.precision(12);
		cout << "generation " << generation << ": the cost of the best solution is: " << best_distance << endl;
	}

	

	//calculate time
	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout.precision(12);
	cout << "the cost of the best solution is: " << best_distance << " in " << totaltime << " seconds" << endl;
	
	//print cost
	currentbest.precision(12);
	currentbest << "cost "<< best_distance << endl;
	//print truck route
	best_chromosome = insert_depot(best_chromosome);
	currentbest<< "1"<<"->";
	for(int l=1;l < best_chromosome.size()-1;l++){
		if(best_chromosome[l] == 0){
			currentbest<< "1"<<endl;
			currentbest<< "1"<<"->";
		}
		else
			currentbest<<best_chromosome[l]+1<<"->";
	}
	currentbest<< "1";

	//print truck route in files with coordinates
	int zeros=0;
	for(int m=1;m < best_chromosome.size();m++){
		if(best_chromosome[m] == 0)
			zeros = zeros + 1;
		if( zeros > 30) break;
	}
	//zeros is number of truck used
	int l=0;
	
	for (int z = 1; z < zeros + 1; z++) {
			string routename = "truck" + num2str(z) + ".txt";
			ofstream routegraph;
			routegraph.open(routename.c_str(), ios::trunc);
			routegraph << individuals[0].num << " " << individuals[0].x << " " << individuals[0].y << endl;
			while (best_chromosome[l] != 0) {
				routegraph << individuals[best_chromosome[l]].num << " " << individuals[best_chromosome[l]].x << " " << individuals[best_chromosome[l]].y << endl;
				l = l + 1;
			}

			routegraph << "1 " << individuals[0].x << " " << individuals[0].y << endl;
			l = l + 1;
			if(l == best_chromosome.size()) break;
		
	}
	

	
	







	return 0;
}



