import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.awt.Color;

public class Fisher {
	
	private int N = 50;
	private double R = 10.;
	private final double D = 1., alpha = 1.;
	private final double dx = 1., dt = .0001;
	
	private double[][] fi = new double[N][N], fiNew = new double[N][N];
	private int[] indexPlus = new int[N], indexMinus = new int[N];
	
	private int visDepth = 5;
	private Visualization vis = new Visualization(N, N);
	
	public Fisher(int N, double R) {
		this.N = N;
		this.R = R;
		this.init();
	}
	
	private void setPlusMinus(){
		indexPlus = new int[N];
		indexMinus = new int[N];
		for(int i=1; i<N; i++){
			indexPlus[i]=i+1;
			indexMinus[i]=i-1;
		}
		indexPlus[N-1] = 0;
		indexMinus[0] = N-1;
	}
	
	private void init(){
		
		this.setPlusMinus();
		fi = new double[N][N];
		double center = N/2.;
		
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++){
				if (Math.sqrt((i-center)*(i-center)+(j-center)*(j-center))<R)
					fi[i][j] = 1.;
				else
					fi[i][j] = 0.;
			}
	}
	
	private void visualize(double[][] grid) {
		for (int col = 0; col < N; col++) 
			for (int row = 0; row < N; row++) 
				vis.set(col, row, grid[col][row] == 0 ? Color.BLACK : Color.getHSBColor((float) (0.666666666666667f*(1f - Math.min(grid[col][row], visDepth)/(float)visDepth)), 1f, 1f));

		vis.draw();
	}
	
	private void setNew(){
		fiNew = new double[N][N];
		
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				fiNew[i][j] = fi[i][j] + D*dt*(fi[i][indexMinus[j]] + 
						fi[i][indexPlus[j]]+fi[indexMinus[i]][j] + 
						fi[indexPlus[i]][j]+fi[i][j])/(dx*dx) + 
						alpha*dt*fi[i][j]*(1-fi[i][j]);
		fi = fiNew;
	}
	
	public void evolve(){
		
		int T = 10000000;
		for(int t=0; t<T; t++){
			this.setNew();
			this.visualize(fi);
		}
	}
	
	/*
	private void setNewJac(){
		
		phiNew = new double[N+1][N+1][N+1];
		
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++)
				for(int k=1; k<N; k++){
					phiNew[i][j][k] = (phi[iPlus[i]][j][k]+phi[iMinus[i]][j][k]+
							phi[i][iPlus[j]][k]+phi[i][iMinus[j]][k]+
							phi[i][j][iPlus[k]]+phi[i][j][iMinus[k]]+
							rho[i][j][k])/6.0;
				}
		
		phi = phiNew;
	}
	
	
	private void setNew(double omega){
		
		//phiNew = new double[N+1][N+1][N+1];
		double phi_ijk_n =0.;
		
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++)
				for(int k=1; k<N; k++){
					phi_ijk_n = phi[i][j][k];
					phi[i][j][k] = (phi[iPlus[i]][j][k]+phi[iMinus[i]][j][k]+
							phi[i][iPlus[j]][k]+phi[i][iMinus[j]][k]+
							phi[i][j][iPlus[k]]+phi[i][j][iMinus[k]]+
							rho[i][j][k])/6.0;
					phi[i][j][k] = (1.-omega)*phi_ijk_n + omega * phi[i][j][k];
				}
		
		//phi = phiNew;
	}
	
	public int evolve(String phiFile, double omega, boolean write) 
			throws FileNotFoundException, UnsupportedEncodingException{
		
		double phi_avg_temp = 0.;
		for(int i=0; i<phi.length; i++)
			for(int j=0; j<phi.length; j++)
				for(int k=0; k<phi.length; k++){
				phi_avg_temp += phi[i][j][k] / 
						(double)(N*N*N);
			}
		
		int count = 0;
		while(true){
			
			this.setNew(omega);
			count++;
			if(count%100==0 && count > 1000){
				double phi_avg = 0.;
				for(int i=1; i<N; i++)
					for (int j=1; j<N; j++)
						for(int k=1; k<N; k++)
							phi_avg += phi[i][j][k] / 
								(double)(N*N*N);
				if(Math.abs((phi_avg_temp - phi_avg)/phi_avg) < 1e-3)
					break;
				else
					phi_avg_temp = phi_avg;
			}
		}
		
		if(write){
			PrintWriter writer1 = new PrintWriter(phiFile, "UTF-8");
			for(int i=1; i<N; i++){
				for(int j=1; j<N; j++)
					for(int k=1; k<N; k++){
						this.getE_ijk(i, j, k);
						writer1.println(i + " " + j + " " + k +  " " + 
									rho[i][j][k] + " " + phi[i][j][k] + " " +
									E_ijk[0] + " " + E_ijk[1] + " " + E_ijk[2]);
					}
				writer1.println();
			}
			writer1.close();
		}
		return count;
		
	}
	
	public void SOR(String SORfile) 
			throws FileNotFoundException, UnsupportedEncodingException{
		
		PrintWriter writer1 = new PrintWriter(SORfile, "UTF-8");
		int Nomega = 50;
		double omega, d_omega = 1./(double)(Nomega);
		int nt = 0;
		for(int w=0; w<Nomega; w++){
			omega = 1. + d_omega * w;
			nt = this.evolve("", omega, false);
			writer1.println(omega + " " + nt);
			this.init();
		}
		writer1.close();
		
	}*/
	
	public static void main(String[] args) 
			throws FileNotFoundException, UnsupportedEncodingException {
		System.out.println("HOLA");
		Fisher fish = new Fisher(50, 10.);
		fish.evolve();
	}
	
}
