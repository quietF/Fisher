import java.awt.Color;

class Walk {
	public static void main(final String[] args) {
		if (args.length != 4) {
			System.err.println("Arguments: width height depth steps/draw");
			return;
		}
		int W = Math.max(1, Integer.parseInt(args[0]));
		int H = Math.max(1, Integer.parseInt(args[1]));
		int D = Math.max(1, Integer.parseInt(args[2]));
		int S = Math.max(1, Integer.parseInt(args[3]));

		int[][] grid = new int[W][H];
		Visualization vis = new Visualization(W, H);

		int col = (int)(0.5*W);
		int row = (int)(0.5*H);
		for (int i = 0; true; i++) {
			grid[col][row]++;
			if (i%S == 0) visualize(vis, grid, W, H, D);

			col = (col + move() + W)%W;
			row = (row + move() + H)%H;
		}
	}

	static void visualize(Visualization vis, int[][] grid, int W, int H, int D) {
		for (int col = 0; col < W; col++) for (int row = 0; row < H; row++) vis.set(col, row, grid[col][row] == 0 ? Color.BLACK : Color.getHSBColor(0.666666666666667f*(1f - Math.min(grid[col][row], D)/(float)D), 1f, 1f));

		vis.draw();
	}

	static int move() {
		double ran = Math.random();
		if (ran < 0.333333333333333) return -1;
		if (ran < 0.666666666666667) return 0;
		else return 1;
	}
}