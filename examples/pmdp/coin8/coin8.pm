//Randomised Consensus Protocol

mdp
const double p1; // in [0.2 , 0.8]
const double p2; // in [0.2 , 0.8]
const double p3; // in [0.2 , 0.8]
const double p4; // in [0.2 , 0.8]
const double p5; 
const double p6;
const double p7;
const double p8;


const int N=8;
const int K;
const int range = 2*(K+1)*N;
const int counter_init = (K+1)*N;
const int left = N;
const int right = 2*(K+1)*N - N;

// shared coin
global counter : [0..range] init counter_init;

module process1
	
	// program counter
	pc1 : [0..3];
	// 0 - flip
	// 1 - write 
	// 2 - check
	// 3 - finished
	
	// local coin
	coin1 : [0..1];	

	// flip coin
	[] (pc1=0)  -> p1 : (coin1'=0) & (pc1'=1) + 1 - p1 : (coin1'=1) & (pc1'=1);
	// write tails -1  (reset coin to add regularity)
	[] (pc1=1) & (coin1=0) & (counter>0) -> (counter'=counter-1) & (pc1'=2) & (coin1'=0);
	// write heads +1 (reset coin to add regularity)
	[] (pc1=1) & (coin1=1) & (counter<range) -> (counter'=counter+1) & (pc1'=2) & (coin1'=0);
	// check
	// decide tails
	[] (pc1=2) & (counter<=left) -> (pc1'=3) & (coin1'=0);
	// decide heads
	[] (pc1=2) & (counter>=right) -> (pc1'=3) & (coin1'=1);
	// flip again
	[] (pc1=2) & (counter>left) & (counter<right) -> (pc1'=0);
	// loop (all loop together when done)
	[done] (pc1=3) -> (pc1'=3);

endmodule

module process2 = process1[pc1=pc2,coin1=coin2,p1=p2] endmodule
module process3 = process1[pc1=pc3,coin1=coin3,p1=p3] endmodule
module process4 = process1[pc1=pc4,coin1=coin4,p1=p4] endmodule
module process5 = process1[pc1=pc5,coin1=coin5,p1=p5] endmodule
module process6 = process1[pc1=pc6,coin1=coin6,p1=p6] endmodule
module process7 = process1[pc1=pc7,coin1=coin7,p1=p7] endmodule
module process8 = process1[pc1=pc8,coin1=coin8,p1=p8] endmodule

label "finished" = pc1=3 &pc2=3 &pc3=3 &pc4=3 & pc5=3 & pc6=3 & pc7=3 & pc8=3;
label "all_coins_equal_1" = coin1 = 1 & coin2 = 1 & coin3 = 1 & coin4 = 1 & coin5 = 1 & coin6 = 1 & coin7 = 1 & coin8 = 1;
label "all_coins_equal_0" = coin1 = 0 & coin2 = 0 & coin3 = 0 & coin4 = 0 & coin5 = 0 & coin6 = 0 & coin7 = 0 & coin8 = 0;
label "agree" = coin1=coin2 & coin2=coin3 & coin3 = coin4 & coin4 = coin5 & coin5 = coin6 & coin6 = coin7 & coin7 = coin8;

rewards "steps"
	true : 1;
endrewards



