import java.util.ArrayList;

public class Matrix {
    private final int rowCount; // 矩阵的行的数量
    private final int colCount; // 矩阵的列的数量
    private final double[][] matrixValue; // 矩阵（值）

    public Matrix(int r, int c) {
        this.rowCount = r;
        this.colCount = c;
        this.matrixValue = new double[r][c];
    }

    public Matrix(int r, int c, double value) {
        this(r, c); // 通过调用另一构造函数避免冗余
        for (int i = 0; i < r; i++) {
            java.util.Arrays.fill(this.matrixValue[i], value); // 使用Arrays.fill来快速填充
        }
    }

    public Matrix(double[][] values) {
        this.rowCount = values.length;
        this.colCount = values[0].length;
        this.matrixValue = new double[rowCount][colCount];
        for (int i = 0; i < rowCount; i++) {
            System.arraycopy(values[i], 0, matrixValue[i], 0, colCount); // 使用System.arraycopy优化复制
        }
    }

    public Matrix(ArrayList<ArrayList<Double>> values) {
        this.rowCount = values.size();
        this.colCount = values.get(0).size();
        this.matrixValue = new double[rowCount][colCount];
        for (int i = 0; i < rowCount; i++) {
            for (int j = 0; j < colCount; j++) {
                matrixValue[i][j] = values.get(i).get(j); // 避免强制转换
            }
        }
    }

    public static void setRow(int r, Matrix row, Matrix matrix) throws Exception {
        if (row.getRowCount() != 1 || row.getColCount() != matrix.getColCount()) {
            throw new Exception("行向量不匹配");
        }
        if (r < 0 || r >= matrix.rowCount) {
            throw new Exception("行数越界");
        }
        for (int i = 0; i < matrix.colCount; i++) {
            matrix.setMatrixValue(r, i, row.getMatrixValue(0, i));
        }
    }

    public static Matrix getEVector(int dimension, int index) {
        Matrix result = new Matrix(1, dimension);
        try {
            result.setMatrixValue(0, index, 1.0);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        return result;
    }

    public static Matrix times(Matrix m1, Matrix m2) throws Exception {
        if (m1.getColCount() != m2.getRowCount()) {
            throw new Exception("矩阵不匹配，无法相乘");
        }

        Matrix result = new Matrix(m1.getRowCount(), m2.getColCount());
        for (int i = 0; i < result.getRowCount(); i++) {
            for (int j = 0; j < result.getColCount(); j++) {
                double sum = 0; // 减少对方法的重复调用
                for (int k = 0; k < m1.getColCount(); k++) {
                    sum += m1.getMatrixValue(i, k) * m2.getMatrixValue(k, j);
                }
                result.setMatrixValue(i, j, sum);
            }
        }
        return result;
    }

    public static Matrix plus(Matrix m1, Matrix m2) throws Exception {
        if (m1.getRowCount() != m2.getRowCount() || m1.getColCount() != m2.getColCount()) {
            throw new Exception("矩阵不匹配，无法相加");
        }

        Matrix result = new Matrix(m1.getRowCount(), m1.getColCount());
        for (int i = 0; i < result.getRowCount(); i++) {
            for (int j = 0; j < result.getColCount(); j++) {
                result.setMatrixValue(i, j, m1.getMatrixValue(i, j) + m2.getMatrixValue(i, j));
            }
        }
        return result;
    }

    public static Matrix dot(Matrix m, double value) throws Exception {
        for (int i = 0; i < m.getRowCount(); i++) {
            for (int j = 0; j < m.getColCount(); j++) {
                double newValue = m.getMatrixValue(i, j) * value;
                m.setMatrixValue(i, j, newValue); // 尽量减少捕获异常的次数
            }
        }
        return m;
    }

    public int getRowCount() {
        return rowCount;
    }

    public int getColCount() {
        return colCount;
    }

    private void checkBounds(int r, int c) throws Exception {
        if (r < 0 || r >= rowCount || c < 0 || c >= colCount) {
            throw new Exception("Out of bound");
        }
    }

    public double getMatrixValue(int r, int c) throws Exception {
        checkBounds(r, c);
        return matrixValue[r][c];
    }

    public void setMatrixValue(int r, int c, double value) throws Exception {
        checkBounds(r, c);
        matrixValue[r][c] = value;
    }

    public Matrix getRow(int r) throws Exception {
        checkBounds(r, 0); // 检查行的有效性
        return new Matrix(new double[][]{matrixValue[r]}); // 使用二维数组构造新的行向量
    }

    public Matrix transpose() throws Exception {
        Matrix transposedMatrix = new Matrix(colCount, rowCount);
        for (int i = 0; i < rowCount; i++) {
            for (int j = 0; j < colCount; j++) {
                transposedMatrix.setMatrixValue(j, i, this.getMatrixValue(i, j));
            }
        }
        return transposedMatrix;
    }
}
