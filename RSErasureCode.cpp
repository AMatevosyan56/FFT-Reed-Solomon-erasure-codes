/*
Encoding/erasure decoding for Reed-Solomon codes over binary extension fields
Author: Sian-Jheng Lin (King Abdullah University of Science and Technology (KAUST), email: sianjheng.lin@kaust.edu.sa)

This program is the implementation of
Lin, Han and Chung, "Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes," FOCS14.
(http://arxiv.org/abs/1404.3458)
*/

#include <cstdlib>
#include <string>
#include <ctime>
#include <cstdint>
#include <cstdio>
#include <iostream>


using GFSymbol = unsigned char;
constexpr uint8_t sizeOfField = 8; // 2^sizeOfField: the size of Galois field
// reducing polynomial for multiplication
GFSymbol mask = 0x1D; // GF(2^8): x^8 + x^4 + x^3 + x^2 + 1  (11101)
GFSymbol Base[sizeOfField] = {1, 214, 152, 146, 86, 200, 88, 230}; // Cantor basis // TODO_ANNA: ???

/*
typedef unsigned short GFSymbol;
#define sizeOfField 16
// reducing polynomial for multiplication
GFSymbol mask = 0x2D;//x^16 + x^5 + x^3 + x^2 + 1 (101101)
GFSymbol Base[sizeOfField] = {1, 44234, 15374, 5694, 50562, 60718, 37196, 16402, 27800, 4312, 27250, 47360, 64952, 64308, 65336, 39198};//Cantor basis
*/

constexpr size_t size = 1<<sizeOfField;
constexpr size_t mod = size - 1;

GFSymbol logTable[size]; // TODO_ANNA: ???
GFSymbol expTable[size]; // TODO_ANNA: ???

//-----Used in decoding procedure-------
GFSymbol skewVec[mod]; // twisted factors used in FFT
GFSymbol B[size>>1]; // factors used in formal derivative
GFSymbol log_walsh[size]; // factors used in the evaluation of the error locator polynomial


//return a*expTable[b] over GF(2^r)
GFSymbol mulE(GFSymbol a, GFSymbol b){
	return a? expTable[(logTable[a]+b&mod) + ((logTable[a]+b) >>sizeOfField)]: 0;
}

void walsh(GFSymbol* data, int size){//fast Walsh–Hadamard transform over modulo mod
	for (int depart_no=1; depart_no<size; depart_no <<= 1){
		for (int j = 0; j < size; j += depart_no<<1){
			for (int i=j; i<depart_no+j; i++){
				unsigned tmp2 = data[i] + mod - data[i+depart_no];
				data[i] = (data[i] + data[i+depart_no]&mod) + (data[i] + data[i+depart_no]>>sizeOfField);
				data[i+depart_no] = (tmp2&mod) + (tmp2>>sizeOfField);
			}
		}
	}
	return;
}

void formal_derivative(GFSymbol* cos, int size){//formal derivative of polynomial in the new basis
	for(int i=1; i<size; i++){
		int leng = ((i^i-1)+1)>>1;
		for(int j=i-leng; j<i; j++)
			cos[j] ^= cos[j+leng];
	}
	for(int i=size; i<size; i<<=1)
		for(int j=0; j<size; j++)
			cos[j] ^= cos[j+i];
	return;
}

// IFFT in the proposed basis
void IFLT(GFSymbol* data, int size, int index){
    // for all layers from the highest to lowest
	for (int depart_no=1; depart_no<size; depart_no <<= 1){
		for (int j=depart_no; j < size; j += depart_no<<1){
			for (int i=j-depart_no; i<j; i++)
                // calculate second set
				data[i+depart_no] ^= data[i];

			GFSymbol skew = skewVec[j+index-1];
			if (skew != mod)
				for (int i=j-depart_no; i<j; i++)
                    // calculate first set
					data[i] ^= mulE(data[i+depart_no], skew);
		}
	}
	return;
}

// FFT in the proposed basis
// index are l shifts in algorithm
void FLT(GFSymbol* data, int size, int index){
    //layers from the lowest to highest
	for(int depart_no = size>>1; depart_no > 0; depart_no >>= 1){
		for (int j = depart_no; j < size; j += depart_no<<1){
            // std::cout<<"depart_no :"<<depart_no<<" ,j :"<<j<<std::endl;
			GFSymbol skew = skewVec[j+index-1];
			if (skew != mod)
				for (int i=j-depart_no; i<j; i++)
                    // calculate first set
					data[i] ^= mulE(data[i+depart_no], skew);
			for (int i=j-depart_no; i<j; i++)
                // calculate second set
				data[i+depart_no] ^= data[i];
		}
        // std::cout<<std::endl;
	}
	return;
}

//initialize logTable[], expTable[]
void init(){
	GFSymbol mas = (1 << (sizeOfField-1))-1;
	GFSymbol state = 1;
	for(int i = 0; i < mod; ++i){
		expTable[state]=i;
         printf("%02X ",state);
        if(state>>(sizeOfField-1)){
        	state &= mas;
                   printf("state now is %02X \n",state);

        	state = state<<1^mask;
                   printf("and now state now is %02X \n",state);

        }else
        	state <<= 1;
    }
    expTable[0] = mod;
    // std::cout<<std::endl;
    // std::cout<<std::endl;
    // std::cout<<std::endl;
    // std::cout<<std::endl;

	for(int i = 0; i < size; ++i){
        printf("%02X ", expTable[i]);

    }
    // std::cout<<std::endl;
    // std::cout<<std::endl;

    logTable[0] = 0;
	for(int i=0; i<sizeOfField; i++)
		for(int j=0; j<1<<i; j++)
			logTable[j+(1<<i)] = logTable[j] ^ Base[i];
    for(int i=0; i<size; i++)
        logTable[i]=expTable[logTable[i]];

    for(int i=0; i<size; i++)
        expTable[logTable[i]]=i;
    expTable[mod] = expTable[0];
    for(int i = 0; i < size; ++i){
               		printf("%02X ", expTable[i]);

    }
    // std::cout<<std::endl;
}


void init_dec(){//initialize skewVec[], B[], log_walsh[]
	GFSymbol base[sizeOfField-1];

	for(int i=1; i<sizeOfField; ++i)
		base[i-1] = 1<<i;

	for(int m=0; m<sizeOfField-1; ++m){
		int step = 1<<(m+1);
		skewVec[(1<<m)-1] = 0;
		for(int i=m; i<sizeOfField-1; i++){
			int s = 1<<(i+1);
			for(int j=(1<<m)-1; j<s; j+=step)
				skewVec[j+s] = skewVec[j] ^ base[i];
		}
		base[m] = mod-logTable[mulE(base[m], logTable[base[m]^1])];
		for(int i=m+1; i<sizeOfField-1; i++)
			base[i] = mulE(base[i], (logTable[base[i]^1]+base[m])%mod);
	}
	for(int i=0; i<size; i++)
		skewVec[i] = logTable[skewVec[i]];

	base[0] = mod-base[0];
	for(int i=1; i<sizeOfField-1; i++)
		base[i] = (mod-base[i]+base[i-1])%mod;

	B[0] = 0;
	for(int i=0; i<sizeOfField-1; i++){
		int depart = 1<<i;
		for(int j=0; j<depart; j++)
			B[j+depart] = (B[j] + base[i])%mod;
	}

	memcpy(log_walsh, logTable, size*sizeof(GFSymbol));
	log_walsh[0] = 0;
	walsh(log_walsh, size);
}

//Encoding algorithm for k/n>0.5: parity is a power of two.
//data: message array. parity: parity array. mem: buffer(size>= n-k)
void encodeH(GFSymbol* knownData, int k, GFSymbol* parity, GFSymbol* mem){
	int t = size-k;
    // TODO_ANNA: check index changes
    // get known k values
    memcpy(mem, knownData, sizeof(GFSymbol)*k);
    IFLT(mem, k, t);
    // keep calculated coefficiets in parity
    for(int j=0; j<k; j++){
        parity[j] = mem[j];
    }
	FLT(parity, k, 0);
	return;
}

//Compute the evaluations of the error locator polynomial
void decode_init(bool* erasure, GFSymbol* log_walsh2){
	for(int i=0; i<size; i++)
		log_walsh2[i] = erasure[i];
	walsh(log_walsh2, size);
	for (int i=0; i<size; i++)
		log_walsh2[i] = (unsigned long)log_walsh2[i]*log_walsh[i]%mod;
	walsh(log_walsh2,size);
	for (int i=0; i<size; i++)
		if(erasure[i]) log_walsh2[i] = mod-log_walsh2[i];
}

void decode_main(GFSymbol* codeword, bool* erasure, GFSymbol* log_walsh2){
	int k2 = size;//k2 can be replaced with k
	for (int i=0; i<size; i++)
        // get F^ = F*P
		codeword[i] = erasure[i]? 0 : mulE(codeword[i], log_walsh2[i]);
	// get coefficients of F^
    IFLT(codeword, size, 0);

    //formal derivative
	for(int i=0; i<size; i+=2){
		codeword[i] = mulE(codeword[i], mod-B[i>>1]);
		codeword[i+1] = mulE(codeword[i+1], mod-B[i>>1]);
	}
	formal_derivative(codeword, k2);
	for(int i=0; i<k2; i+=2){
		codeword[i] = mulE(codeword[i], B[i>>1]);
		codeword[i+1] = mulE(codeword[i+1], B[i>>1]);
	}

    // get values of format derivative
	FLT(codeword, k2, 0);
	for (int i=0; i<k2; i++)
		codeword[i] = erasure[i]? mulE(codeword[i], log_walsh2[i]) : 0;
}


//TODO_ANNA: separate into functions all preparation part.
void test(int k){
	//-----------Generating message----------
	GFSymbol data[size] = {0}; // message array
	srand(time(NULL));
	for(int i=size-k; i<size; ++i)
		data[i] = rand()&mod; //fill with random numbers

	printf("Message(First n-k are zeros): \n");
	for(int i=0; i<size; ++i)
		printf("%02X ", data[i]);
	printf("\n");

	//---------encoding----------
	GFSymbol codeword[size];
	encodeH(&data[size-k], k, data, codeword);

	memcpy(codeword, data, sizeof(GFSymbol)*size);

	printf("Codeword:\n");
	for(int i=0; i<size; i++)
		printf("%02X ", codeword[i]);
	printf("\n");

	//--------erasure simulation---------
	bool erasure[size] = {0};//Array indicating erasures
	for(int i=k; i<size; ++i)
		erasure[i] = 1;

    // randomly mix putted 1s in erasue to simplate random erasion of n-k elements
    // permuting the erasure array
	for(int i=size-1; i>0; i--){
		int pos = rand()%(i+1);
		if(i != pos){
			bool tmp = erasure[i];
			erasure[i] = erasure[pos];
			erasure[pos] = tmp;
		}
	}

    // erase selected places of codeword
	for (int i=0; i<size; i++)//erasure codeword symbols
		if(erasure[i]) codeword[i] = 0;

	printf("Erasure (XX is erasure):\n");
	for(int i=0; i<size; i++){
		if(erasure[i]) printf("XX ");
		else printf("%02X ", codeword[i]);
	}
	printf("\n");


	//---------Erasure decoding----------------
	GFSymbol log_walsh2[size];
    // in log_walsh2 there are values of error locating polynomial over field
	decode_init(erasure, log_walsh2);//Evaluate error locator polynomial
	//---------main processing----------
	decode_main(codeword, erasure, log_walsh2);

	printf("Decoded result:\n");
	for(int i=0; i<size; i++){
		if(erasure[i]) printf("%02X ", codeword[i]);
		else printf("XX ");
	}
	printf("\n");

	for (int i=0; i<size; i++){//Check the correctness of the result
		if(erasure[i] == 1)
			if(data[i] != codeword[i]){
				printf("Decoding Error!\n");
				return;
			}
	}
	printf("Decoding is successful!\n");
	return;
}

int main(){
	init(); //fill log table and exp table
	init_dec(); //compute factors used in erasure decoder
	test(size/2); //test(int k), k: message size
	return 0;
}
