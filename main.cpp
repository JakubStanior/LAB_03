#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <complex>
#include <iomanip>
#include <algorithm>

using namespace std;

const int A = 7;
const int B = 0;
const int C = 5;
const int D = 4;
const int E = 4;
//const int F;
const int ARRSIZE = 715; // ilosc probek: <0; ABC>

const double pi = std::acos(-1);
const complex<double> im(0, 1);

double editSampleOfSignal(double t, char funcType)
{
	int fn=4 ;// czestotliwosc dla sygnalu s
	double sum;
	float amp=1.0;
	float fi=C*M_PI;
	switch (funcType)
	{
		case 's':
			return amp*sin(2 * M_PI * fn * t + fi);
		break;
		
		case 'x':
			return A*t*t + B*t + C;
		break;
		
		case 'y':
			return 2 * (editSampleOfSignal(t, 'x')*editSampleOfSignal(t, 'x')) + 12*cos(t);
		break;
		
		case 'z':
			return sin(2*M_PI*7*t) * editSampleOfSignal(t, 'x') - 0.2 * log10(abs(editSampleOfSignal(t, 'y')) + M_PI);
		break;
		
		case 'u':
			return sqrt(abs(editSampleOfSignal(t, 'y')*editSampleOfSignal(t, 'y')*editSampleOfSignal(t, 'z') )) - 1.8 * sin(0.4 * t * editSampleOfSignal(t, 'z') * editSampleOfSignal(t, 'x'));
		break;
		
		case 'v':
			if(0.22 > t && t >= 0) return (1-7*t) * sin( (2*M_PI*t*10) / (t+0.04) );
			if(0.7 > t && t >= 0.22) return 0.63 * t * sin(125*t);
			if(1.0 >= t && t >= 0.7) return pow(t,-0.662) + 0.77 * sin(8*t);
			if(1.0 < t) return 0;
		break;
		
		case '1':
			sum=0;
			for(int n=1;n<=2;n++)
			{
				sum+= (cos(12*t*n*n) + cos(16*t * n)) / (n*n);
			}
			return sum;
		break;
		
		case '2':
			sum=0;
			for(int n=1;n<=4;n++)
			{
				sum+= (cos(12*t*n*n) + cos(16*t * n)) / (n*n);
			}
			return sum;
		break;
		
		case '3':
			sum=0;
			//n<= A B
			for(int n=1;n<=71;n++)
			{
				sum+= (cos(12*t*n*n) + cos(16*t * n)) / (n*n);
			}
			return sum;
		break;
		
		default:
			cout << "Brak takiej funkcji" << endl;
		break;
	}
}

double* generateSignal(double* time_stamps, char funcType)
{
	double* result = new double[ARRSIZE];
	for (int i=0;i<ARRSIZE;i++)
	{
		double tmp = editSampleOfSignal(time_stamps[i], funcType);
		//if (abs(tmp)<pow(10,-5))result[i] = 0;
		result[i]=tmp;
	}
	return result;
}

complex<double>* dft(double* f)
{
	complex<double>* result = new complex<double>[ARRSIZE];
	for(int k=0;k<ARRSIZE;++k)
	{
		complex<double> sum=0;
		for(int n=0;n<ARRSIZE;n++)
		{
			sum+=f[n]*exp(-2.0* M_PI * im *(double)k*(((double)n)/ARRSIZE));
		}
		result[k]=sum;
	}
	return result;
}

double* spectrum(complex<double>* f)
{
	double* result = new double[ARRSIZE];
	for(int i=0; i<ARRSIZE; i++)
	{
		result[i] = sqrt((imag(f[i])*imag(f[i])) + (real(f[i])*real(f[i])));
		//result[i] = abs(f[i]);
	}
	
	return result;
}

void zeroThreshold(double* func)
{
	double max=func[0];
	for(int i=0; i<ARRSIZE; i++)
	{
		if(func[i] > max) max = func[i];
	}
	
	double threshold = max / 10000;
	
	for(int i=0; i<ARRSIZE; i++)
	{
		if(func[i] < threshold) func[i]=0;
	}
}

double* convertToDecibelScale(double* func)
{
	double* result = new double[ARRSIZE];
	for(int i=0; i<ARRSIZE; i++)
	{
		result[i] = 10.0 * log10(func[i]) + pow(10,(-6));
		if(result[i]<0) result[i]=0;
	}
	return result;
}

double* idft(complex<double>* ComplexArray)
{
	double tmp;
	double* result = new double[ARRSIZE];
	for(int n=0; n<ARRSIZE; n++)
	{
		tmp=0;
		for(int k=0;k<ARRSIZE;k++)
		{
			double phase = (2 * M_PI / ARRSIZE) * k * n;
			tmp += cos(phase) * ComplexArray[k].real() - sin(phase) * ComplexArray[k].imag();
		}
		tmp /= ARRSIZE;
		result[n]=tmp;
	}
	return result;
}

void writeToFile(double* x, double* y, char func, int numberOfPlot)
{
	string fileName(1, func);
	if(fileName=="1" || fileName=="2" || fileName=="3") fileName= "p" + fileName;
	
	fileName += "-";
	fileName += to_string(numberOfPlot);
	fileName += ".csv";
	
	ofstream file1;
	file1.open(fileName.c_str());
	if(file1.good())
	{
		for(int i=0;i<ARRSIZE; i++)
		{
			file1 << x[i] << "," << y[i] << endl;
		}
	}
	else cout << "Error opening " << func << ".csv"<<endl;
	file1.close();
	cout << "Generated " << fileName << endl;
}

void generateFileFromSignal(double* time_stamp, char funcType, int fs)
{
	double* signal_data1 = new double[ARRSIZE];
	double* signal_data2 = new double[ARRSIZE];
	double* frequencyScale = new double[ARRSIZE];
	signal_data1 = generateSignal(time_stamp, funcType);
	writeToFile(time_stamp, signal_data1, funcType, 1);
	
	//ETAP 2
	complex<double>* cDFT = new complex<double>[ARRSIZE];
	cDFT = dft(signal_data1);
	signal_data1 = spectrum(cDFT);
	
	for(int i=0; i<ARRSIZE; i++)
	{
		frequencyScale[i]=(i*(double)fs/(double)ARRSIZE);
		if(i<ARRSIZE/2) signal_data1[i] = (signal_data1[i] * 2.0) / ARRSIZE;
		else 
		{
			signal_data1[i]=0;
			frequencyScale[i]=0;
		}
	}
	writeToFile(frequencyScale, signal_data1, funcType, 2);
	
	//ETAP3
	signal_data2 = spectrum(cDFT);
	signal_data2 = convertToDecibelScale(signal_data2);
	
	for(int i=0; i<ARRSIZE; i++)
	{
		if(i>ARRSIZE/2) signal_data2[i]=0;
	}
	
	// wykres3
	writeToFile(frequencyScale, signal_data2, funcType, 3); 
	cout << "----------------------------------------" << endl;
	
	delete [] signal_data1, signal_data2, frequencyScale;
}

void generateIDFT(double *time_stamp)
{
	double* signal_data = new double[ARRSIZE];
	signal_data = generateSignal(time_stamp, 's');
	writeToFile(time_stamp, signal_data, 'f', 1);
	complex<double>* cDFT = new complex<double>[ARRSIZE];
	cDFT = dft(signal_data);
	signal_data = idft(cDFT);
	writeToFile(time_stamp, signal_data, 'r', 1);
}

int main(int argc, char** argv) 
{	
	const int fs = 1000;  // Czestotliwosc probkowania
	const double dt= 1.0/ fs; 	// Krok
	
	double* time_stamp = new double[ARRSIZE];
	for (int i=0;i<ARRSIZE; i++)
	{
		time_stamp[i] = i*dt;
	}
	generateFileFromSignal(time_stamp, 's', fs);
	generateFileFromSignal(time_stamp, 'x', fs);
	generateFileFromSignal(time_stamp, 'y', fs);
	generateFileFromSignal(time_stamp, 'z', fs);
	generateFileFromSignal(time_stamp, 'u', fs);
	generateFileFromSignal(time_stamp, 'v', fs);
	generateFileFromSignal(time_stamp, '1', fs);
	generateFileFromSignal(time_stamp, '2', fs);
	generateFileFromSignal(time_stamp, '3', fs);
	//Zad1, 4
	generateIDFT(time_stamp);
	
	delete [] time_stamp;
	return 0;
}

