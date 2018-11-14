Func void main()
{
	Matrix data,Data,data2,Data2;
	Matrix leg1,leg2,leg3,leg4;
	Matrix body;
	Real L1,L2,L3;
	Integer i,j,k;
	Integer col,d_num,d_col;
	Integer win;
	
	//read data << "3.mat";
	read data << "logging.mat";
	Data = trans(data);
	//read data2 << "2-3.mat";
	//Data2 = trans(data2);
	
	L1 = 0.130;	L2 = 0.03677;	L3 = 0.061;
	
	col = Cols(Data);	d_num = 6;	d_col = col/d_num;
	
	for(i=0;i<=d_num;i++)
	{
		leg1(i*2+1,1) = Data(2+3,i*d_col+1);
		leg1(i*2+1,2)= leg1(i*2+1,1)-L3*cos(Data(6+3,i*d_col+1));
		leg1(i*2+1,3)= leg1(i*2+1,2)-L2*cos(Data(5+3,i*d_col+1));
		leg1(i*2+1,4)= leg1(i*2+1,3)-L1*cos(Data(4+3,i*d_col+1));
		leg1(i*2+2,1)= Data(3+3,i*d_col+1);
		leg1(i*2+2,2)= leg1(i*2+2,1)-L3*sin(Data(6+3,i*d_col+1));
		leg1(i*2+2,3)= leg1(i*2+2,2)-L2*sin(Data(5+3,i*d_col+1));
		leg1(i*2+2,4)= leg1(i*2+2,3)-L1*sin(Data(4+3,i*d_col+1));
	}
	for(i=0;i<=d_num;i++)
	{
		leg2(i*2+1,1) = Data(7+3,i*d_col+1);
		leg2(i*2+1,2)= leg2(i*2+1,1)-L3*cos(Data(11+3,i*d_col+1));
		leg2(i*2+1,3)= leg2(i*2+1,2)-L2*cos(Data(10+3,i*d_col+1));
		leg2(i*2+1,4)= leg2(i*2+1,3)-L1*cos(Data(9+3,i*d_col+1));
		leg2(i*2+2,1)= Data(8+3,i*d_col+1);
		leg2(i*2+2,2)= leg2(i*2+2,1)-L3*sin(Data(11+3,i*d_col+1));
		leg2(i*2+2,3)= leg2(i*2+2,2)-L2*sin(Data(10+3,i*d_col+1));
		leg2(i*2+2,4)= leg2(i*2+2,3)-L1*sin(Data(9+3,i*d_col+1));
	}
	for(i=0;i<=d_num;i++)
	{
		leg3(i*2+1,1) = Data(12+3,i*d_col+1);
		leg3(i*2+1,2)= leg3(i*2+1,1)-L3*cos(Data(16+3,i*d_col+1));
		leg3(i*2+1,3)= leg3(i*2+1,2)-L2*cos(Data(15+3,i*d_col+1));
		leg3(i*2+1,4)= leg3(i*2+1,3)-L1*cos(Data(14+3,i*d_col+1));
		leg3(i*2+2,1)= Data(13+3,i*d_col+1);
		leg3(i*2+2,2)= leg3(i*2+2,1)-L3*sin(Data(16+3,i*d_col+1));
		leg3(i*2+2,3)= leg3(i*2+2,2)-L2*sin(Data(15+3,i*d_col+1));
		leg3(i*2+2,4)= leg3(i*2+2,3)-L1*sin(Data(14+3,i*d_col+1));
	}
	for(i=0;i<=d_num;i++)
	{
		leg4(i*2+1,1) = Data(17+3,i*d_col+1);
		leg4(i*2+1,2)= leg4(i*2+1,1)-L3*cos(Data(21+3,i*d_col+1));
		leg4(i*2+1,3)= leg4(i*2+1,2)-L2*cos(Data(20+3,i*d_col+1));
		leg4(i*2+1,4)= leg4(i*2+1,3)-L1*cos(Data(19+3,i*d_col+1));
		leg4(i*2+2,1)= Data(18+3,i*d_col+1);
		leg4(i*2+2,2)= leg4(i*2+2,1)-L3*sin(Data(21+3,i*d_col+1));
		leg4(i*2+2,3)= leg4(i*2+2,2)-L2*sin(Data(20+3,i*d_col+1));
		leg4(i*2+2,4)= leg4(i*2+2,3)-L1*sin(Data(19+3,i*d_col+1));
	}
	for(i=0;i<=d_num;i++)
	{
		body(i*2+1,1)	= leg1(i*2+1,4);
		body(i*2+1,2)	= leg2(i*2+1,4);
		body(i*2+1,3)	= leg3(i*2+1,4);
		body(i*2+1,4)	= leg4(i*2+1,4);
		body(i*2+1,5)	= leg1(i*2+1,4);
		body(i*2+2,1)	= leg1(i*2+2,4);
		body(i*2+2,2)	= leg2(i*2+2,4);
		body(i*2+2,3)	= leg3(i*2+2,4);
		body(i*2+2,4)	= leg4(i*2+2,4);
		body(i*2+2,5)	= leg1(i*2+2,4);
	}
	
	win = 1;
	
	///*
	mgplot_cmd(win,"set zeroaxis");
	mgplot_cmd(win,"set xrange [0.0:7]");
	mgplot_cmd(win,"set yrange [0.5:1.5]");
	//mgplot_cmd(win,"set xrange [0:7]");
	//mgplot_cmd(win,"set yrange [0:2]");
	mgplot_cmd(win,"set xtics 1.0");
	mgplot_cmd(win,"set ytics 0.5");
	mgplot_cmd(win,"set grid");
	mgplot_cmd(win,"set ylabel'Y [m]'");
	mgplot_cmd(win,"set xlabel'X [m]'");
	mgplot_cmd(win,"set border 31 lw 2");
	mgplot_cmd(win,"set key below");
	mgplot_cmd(win,"set size 1.5,0.7");
	mgplot(win,Data(2+3,:),Data(3+3,:),{""},{"lw 6"});
	mgreplot(win,Data(7+3,:),Data(8+3,:),{""},{"lw 6"});
	mgreplot(win,Data(12+3,:),Data(13+3,:),{""},{"lw 6"});
	mgreplot(win,Data(17+3,:),Data(18+3,:),{""},{"lw 6"});
	pause(1.0);
	for(i=0;i<=d_num;i++)
	{
		mgreplot(win,body(i*2+1,:),body(i*2+2,:),{""},{"lt 7 lw 4"});
		pause(1.0);
		mgreplot(win,leg1(i*2+1,:),leg1(i*2+2,:),{""},{"lt 7 lw 4"});
		pause(1.0);
		mgreplot(win,leg2(i*2+1,:),leg2(i*2+2,:),{""},{"lt 7 lw 4"});
		pause(1.0);
		mgreplot(win,leg3(i*2+1,:),leg3(i*2+2,:),{""},{"lt 7 lw 4"});
		pause(1.0);
		mgreplot(win,leg4(i*2+1,:),leg4(i*2+2,:),{""},{"lt 7 lw 4"});
		pause(1.0);
	}
	mgplot_cmd(win,"set terminal post eps color solid enhanced 25");
	mgplot_cmd(win,"set output 'x-y.eps'");
	mgplot_replot(win);
	win++;
	//*/
	
	/*
	mgplot_cmd(win,"set zeroaxis");
	mgplot_cmd(win,"set grid");
	mgplot_cmd(win,"set ylabel'Velocity [m/s]'");
	mgplot_cmd(win,"set xlabel'Time [s]'");
	mgplot_cmd(win,"set key below");
	mgplot(win,Data(1,:),Data(22,:),{"leg1"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(26,:),{"leg2"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(30,:),{"leg3"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(34,:),{"leg4"},{"lw 6"});
	mgplot_cmd(win,"set terminal post eps color solid enhanced 25");
	mgplot_cmd(win,"set output 'Velocity.eps'");
	mgplot_replot(win);
	win++;
	//*/
	
	/*
	mgplot_cmd(win,"set zeroaxis");
	mgplot_cmd(win,"set grid");
	mgplot_cmd(win,"set ylabel'Angular velocity [rad/s]'");
	mgplot_cmd(win,"set xlabel'Time [s]'");
	mgplot_cmd(win,"set key below");
	mgplot(win,Data(1,:),Data(23,:),{"leg1-1"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(24,:),{"leg1-2"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(25,:),{"leg1-3"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(27,:),{"leg2-1"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(28,:),{"leg2-2"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(29,:),{"leg2-3"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(31,:),{"leg3-1"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(32,:),{"leg3-2"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(33,:),{"leg3-3"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(35,:),{"leg4-1"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(36,:),{"leg4-1"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(37,:),{"leg4-1"},{"lw 6"});
	mgplot_cmd(win,"set terminal post eps color solid enhanced 25");
	mgplot_cmd(win,"set output 'Angular_velocity.eps'");
	mgplot_replot(win);
	win++;
	//*/
	
	/*
	mgplot_cmd(win,"set zeroaxis");
	mgplot_cmd(win,"set grid");
	mgplot_cmd(win,"set ylabel'c_error [N.D.]'");
	mgplot_cmd(win,"set xlabel'Time [s]'");
	mgplot_cmd(win,"set key below");
	mgplot(win,Data(1,:),Data(38,:),{"c_error1"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(39,:),{"c_error2"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(40,:),{"c_error3"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(41,:),{"c_error4"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(42,:),{"c_error5"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(43,:),{"c_error6"},{"lw 6"});
	mgplot_cmd(win,"set terminal post eps color solid enhanced 25");
	mgplot_cmd(win,"set output 'c_error.eps'");
	mgplot_replot(win);
	win++;
	//*/
	
	/*
	mgplot_cmd(win,"set zeroaxis");
	mgplot_cmd(win,"set grid");
	mgplot_cmd(win,"set ylabel'||F||'");
	mgplot_cmd(win,"set xlabel'Time [s]'");
	mgplot_cmd(win,"set key below");
	mgplot(win,Data(1,:),Data(62,:),{"leg1"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(63,:),{"leg2"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(64,:),{"leg3"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(65,:),{"leg4"},{"lw 6"});
	mgplot_cmd(win,"set terminal post eps color solid enhanced 25");
	mgplot_cmd(win,"set output 'cgmres_error.eps'");
	mgplot_replot(win);
	win++;
	//*/
	
	/*
	mgplot_cmd(win,"set zeroaxis");
	mgplot_cmd(win,"set grid");
	mgplot_cmd(win,"set ylabel'optimization time [ms]'");
	mgplot_cmd(win,"set xlabel'Time [s]'");
	mgplot_cmd(win,"set key below");
	mgplot(win,Data(1,:),Data(66,:),{"leg1"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(66,:),{"leg2"},{"lw 6"});
	mgreplot(win,Data(1,:),Data(67,:),{"leg3"},{"lw 6"});
	mgreplot(win,Data2(1,:),Data2(67,:),{"leg4"},{"lw 6"});
	mgplot_cmd(win,"set terminal post eps color solid enhanced 25");
	mgplot_cmd(win,"set output 'optimization time.eps'");
	mgplot_replot(win);
	win++;
	//*/
	
	pause;
}
main()