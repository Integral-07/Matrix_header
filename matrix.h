#pragma once

#include<iostream>
#include<iomanip>
#include<string>


static unsigned int PREC = 3; //有効桁数


template<class Type>
class BaseMatrix {
	//Matrix,Determinantの基底クラス

protected:

	Type** array = nullptr;
	int row = 0, col = 0;

	// すべて0の行があるかチェック
	bool check_AllZero() {

		int bit = 0;

		for (int i = 0; i < row; i++) {

			bit = 0; //符号ビットの初期化
			for (int j = 0; j < col; j++) {

				if (array[i][j]);
				else bit++;
			}

			if (!bit) return true;
		}

		return false;
	}


public:

	//コンストラクタ
	explicit BaseMatrix(Type* x, int _row, int _col) : row(_row), col(_col) {

		//array配列を動的確保
		array = new Type * [row];
		for (int i = 0; i < row; i++) {

			array[i] = new Type[col];
		}

		//配列をコピー
		for (int c = 0; c < row; c++) {
			for (int d = 0; d < col; d++) {

				array[c][d] = x[d + c * col];
			}
		}
	}

	//デフォルトコンストラクタ
	explicit BaseMatrix() : row(2), col(2) {

		Type temp[4] = { 1, 0, 0, 1 };

		//array配列を動的確保
		array = new Type * [row];

for (int i = 0; i < row; i++) {
	array[i] = new Type[col];
}

//配列をコピー
for (int c = 0; c < row; c++) {
	for (int d = 0; d < col; d++) {

		array[c][d] = temp[c * col + d];
	}
}
	}

	//コピーコンストラクタ
	BaseMatrix(const BaseMatrix& bm) {

		if (&bm != this) {
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < col; j++) {

					this->array[i][j] = bm.get_element(i, j);
				}
			}

			this->row = bm.get_row();
			this->col = bm.get_col();
		}
	}

	//仮想デストラクタ
	virtual ~BaseMatrix() {

		for (int i = 0; i < this->row; i++) {

			delete[] this->array[i];
		}
		delete[] this->array;
	}

	//ゲッタ
	int get_row() const { return row; }										//行を取得
	int get_col() const { return col; }										//列を取得
	Type get_element(int _row, int _col) const { return array[_row][_col]; }//要素を取得


	//単位行列にする
	void to_Identity() {

		if (row == col) {

			for (int i = 0; i < row; i++)
				for (int j = 0; j < col; j++) {

					if (i == j) {
						array[i][j] = 1;
					}
					else {
						array[i][j] = 0;
					}
				}
		}

		else {

			std::cout << "型が正しくありません." << std::endl;
		}
	}


	//ゼロ行列にする
	void to_Zero() {

		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++) {

				array[i][j] = 0;
			}
	}


	//掃き出し(オブジェクトに変更を加えます)
	void Sweep() {

		if (check_AllZero()) {

			std::cout << "要素がすべてゼロの行があります。正しく掃き出せていない可能性があります。\n";
		}

		Type K, L;

		for (int i = 0; i < row; i++) {


			//i列目の掃き出し
			K = array[i][i];
			for (int j = i; j < col; j++) {

				if (K == 0 || i == row - 1) break;

				array[i][j] /= K;

			}

			for (int r = i + 1; r < row; r++) {
				L = array[r][i];

				for (int a = 0; a < col; a++) {

					if (L == 0) break;

					array[r][a] += array[i][a] * -L;

				}
			}
		}
	}

};


template<class Type = double>
class Matrix : public BaseMatrix<Type> {
	//行列を扱うクラス

private:

	//配列要素の積の積和を返す
	friend Type Sigma_array(const Matrix<Type>& x, const Matrix<Type>& y, int _row, int _col) {

		if (x.col == y.row) {

			Type result = 0;

			for (int k = 0; k < x.col; k++) {
				result += x.array[_row][k] * y.array[k][_col];
			}

			return result;
		}

		else {
			return -1;
		}
	}

	//add関数
	Matrix<Type> add(const Matrix& x, const Matrix& y) const;

	//sub関数
	Matrix<Type> sub(const Matrix& x, const Matrix& y) const;

	//pro関数
	Matrix<Type> pro(const Matrix& x, const Matrix& y) const;

	//pro_n関数
	Matrix<Type> pro_n(const Matrix& x, Type nd) const;

public:

	//コンストラクタ
	explicit Matrix(Type* x, int _row, int _col) : BaseMatrix<Type>(x, _row, _col) { }

	//デフォルトコンストラクタ
	explicit Matrix() : BaseMatrix<Type>() { }

	//デストラクタ
	~Matrix() { }

	Matrix T(); //転置行列にする

	bool Is_regular(); //正則判定

	Type cof(int _i = 0, int _j = 0); //余因子を求める

	Matrix cof_matrix(); // 余因子行列を返す

	Matrix Inverse(); //逆行列を返す

	int Rank(); //階級を求める

	Matrix Zero(int _row, int _col) const; //ゼロ行列を生成する

	Matrix Identity(int _row, int _col) const; //単位行列を生成する

	

//演算子のオーバーロード

	//==演算子のオーバーロード
	friend bool operator==(const Matrix<Type>& x, const Matrix<Type>& y) {
		//実数の等価比較は誤差を考慮できないので、要素が整数のときの使用のみを想定する

		if (x.row == y.row && x.col == y.col) {

			for (int i = 0; i < x.row; i++) {
				for (int j = 0; j < x.col; j++) {

					if ((int)x.array[i][j] == (int)y.array[i][j])
						continue;

					else
						return false;
				}
			}

			return true;
		}

		else {

			return false;
		}
	}
	
	//!=演算子のオーバーロード
	friend bool operator!=(const Matrix<Type>& x, const Matrix<Type>& y) {
		//実数の等価比較は誤差を考慮できないので、要素が整数のときの使用のみを想定する

		return !(x == y);
	}

	// |演算子のオーバーロード(行列の結合)
	friend Matrix<Type> operator|(const Matrix<Type>& x, const Matrix<Type>& y) {

		if (x.row == y.row && typeid(x) == typeid(y)) {

			Type* temp_array = nullptr;
			temp_array = new Type[x.row * (x.col + y.col)];

			for (int i = 0; i < x.row; i++) {
				for (int s = 0; s < x.col; s++) {

					temp_array[i * (x.col + y.col) + s] = x.array[i][s];
				}

				for (int t = 0; t < y.col; t++) {

					temp_array[i * (x.col + y.col) + t + x.col] = y.array[i][t];
				}
			}

			return Matrix(temp_array, x.row, (x.col + y.col));
		}

		else {

			std::cout << "型が正しくありません" << std::endl;
			return x.Zero(1, 1);
		}
	}

	//+=演算子のオーバーロード
	Matrix operator+=(const Matrix& x) {

		return add(*this, x);
	}

	//-=演算子のオーバーロード
	Matrix operator-=(const Matrix& x) {

		return sub(*this, x);
	}

	//*演算子のオーバーロード
	Matrix operator*(Type nd) {

		return pro_n(*this, nd);
	}

	//*=演算子のオーバーロード
	Matrix operator*=(Type nd) {

		return pro_n(*this, nd);
	}

	//*=演算子のオーバーロード
	Matrix operator*=(const Matrix& x) {

		return pro(*this, x);
	}

	//+演算子のオーバーロード
	friend Matrix operator+(const Matrix& x, const Matrix& y) {

		return x.add(x, y);
	}

	//-演算子のオーバーロード
	friend Matrix operator-(const Matrix& x, const Matrix& y) {

		return x.sub(x, y);
	}

	//*演算子のオーバーロード
	friend Matrix operator*(const Matrix& x, const Matrix& y) {

		return x.pro(x, y);
	}

	//挿入子のオーバーロード
	friend std::ostream& operator<<(std::ostream& s, Matrix& x) {

		int digits;
		int max_digits = 0;
		Type value;
		int prec = 0;

		for (int n = 0; n < x.row; n++) {
			for (int m = 0; m < x.col; m++) {

				value = x.array[n][m];
				digits = 0;

				if (value - (int)value != 0) {
					//小数時の処理

					prec = PREC;   //有効桁数
					digits += prec + 2;
				}

				else {

					std::string d = std::to_string((int)value);


					digits += d.size() / 8 + d.size() % 8;
				}


				if (max_digits < digits) max_digits = digits;
			}
		}

		for (int i = 0; i < x.row; i++) {

			if (i == 0) {

				s << "_" << std::setw((max_digits + 1) * x.col + 3) << "_\n";
			}

			s << "| ";
			for (int j = 0; j < x.col; j++) {

				s.precision(prec);
				s.width(max_digits);
				s << x.array[i][j] << " ";
			}
			s << "|\n";

		}

		s << "-" << std::setw((max_digits + 1) * x.col + 3) << "-\n";

		return s;
	}

};


template<class Type = double>
class Determinant : public BaseMatrix<Type> {
	//行列式を扱うクラス

private:

	int pointed_row, pointed_col;

	//余因子(行列式計算用)
	Determinant<Type> cof_for_det(int i = 0, int j = 0) const; 

	Type add(const Determinant& x, const Determinant& y) const {

		return x.det() + y.det();
	}

	Type sub(const Determinant& x, const Determinant& y) const {

		return x.det() - y.det();
	}

	Type pro(const Determinant& x, const Determinant& y) const {

		return x.det() * y.det();
	}

	Determinant pro_n(const Determinant& x, Type nd) {

		Type* temp_array;
		temp_array = new Type[x.row * x.col];

		for (int i = 0; i < x.row; i++) {
			for (int j = 0; j < x.col; j++) {

				if (pointed_row == i || pointed_col == j) {

					temp_array[i * x.col + j] = (x.array[i][j]) * nd;
				}
				else {

					temp_array[i * x.col + j] = x.array[i][j];
				}
			}
		}

		return	Determinant(temp_array, x.row, x.col);
	}

public:

	//コンストラクタ
	explicit Determinant(Type* x, int _row, int _col) : BaseMatrix<Type>(x, _row, _col), pointed_row(-1), pointed_col(-1){ }

	//デフォルトコンストラクタ
	explicit Determinant() : BaseMatrix<Type>(), pointed_row(-1), pointed_col(-1) { }

	//デストラクタ
	~Determinant() { }


	Determinant<Type> Zero(int _row, int _col) const; //要素がすべてゼロの行列式を生成

	Determinant<Type> Identity(int _row, int _col) const; //要素が単位行列の行列式を生成

	Type det() const; //行列式を求める

	void set_point(int point2row = -1, int point2col = -1); //演算対象の行、列を指定する(実数倍するとき)
	

//演算子のオーバーロード

	//*演算子のオーバーロード(実数倍)
	Determinant operator*(Type nd) {

		return pro_n(*this, nd);
	}

	//*=演算子のオーバーロード(実数倍)
	Determinant operator*=(Type nd) {

		return pro_n(*this, nd);
	}

	//+演算子のオーバーロード
	friend Type operator+(const Determinant& x, const Determinant& y) {

		return x.add(x, y);
	}

	//-演算子のオーバーロード
	friend Type operator-(const Determinant& x, const Determinant& y) {

		return x.sub(x, y);
	}

	//*演算子のオーバーロード
	friend Type operator*(const Determinant& x, const Determinant& y) {

		return x.pro(x, y);
	}

	//==演算子のオーバーロード
	friend bool operator==(const Determinant<Type>& x, const Determinant<Type>& y) {

		if (x.det() == y.det()) 
			return true;

		else 
			return false;
	}

	//!=演算子のオーバーロード
	friend bool operator!=(const Determinant<Type>& x, const Determinant<Type>& y) {

		return !(x == y);
	}

	//挿入子のオーバーロード
	friend std::ostream& operator<<(std::ostream& s, const Determinant& x) {

		int digits;
		int max_digits = 0;
		Type value;
		int prec = 1;

		for (int n = 0; n < x.row; n++) {
			for (int m = 0; m < x.col; m++) {

				digits = 0;
				value = x.array[n][m];

				if (value - (int)value != 0) {
					//小数時の処理

					prec = PREC;   //有効桁数
					digits += prec + 2;
				}

				else {

					std::string d = std::to_string((int)value);


					digits += d.size() / 8 + d.size() % 8;
				}


				if (max_digits < digits) max_digits = digits;
			}
		}

		for (int i = 0; i < x.row; i++) {

			s << "| ";
			for (int j = 0; j < x.col; j++) {

				s.precision(prec);
				s.width(max_digits);
				s << x.array[i][j] << " ";
			}
			s << "|\n";

		}

		return s;
	}

};