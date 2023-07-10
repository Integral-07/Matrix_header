#include "Matrix.h"


template class Matrix<double>;
template class Matrix<int>;


//add関数
template<class Type>
Matrix<Type> Matrix<Type>::add(const Matrix& x, const Matrix& y) const {

	Type* temp_array = nullptr;
	temp_array = new Type[x.row * x.col];

	if (x.row == y.row && x.col == y.col) {

		for (int i = 0; i < x.row; i++) {
			for (int j = 0; j < x.col; j++) {
				temp_array[i * x.col + j] = x.array[i][j] + y.array[i][j];
			}
		}
		return Matrix(temp_array, x.row, x.col);
	}

	else {

		std::cout << "型が正しくありません." << std::endl;
		return x.Zero(1, 1);
	}
}

//sub関数
template<class Type>
Matrix<Type> Matrix<Type>::sub(const Matrix& x, const Matrix& y) const {

	Type* temp_array = nullptr;
	temp_array = new Type[x.row * x.col];

	if (x.row == y.row && x.col == y.col) {

		for (int i = 0; i < x.row; i++) {
			for (int j = 0; j < x.col; j++) {

				temp_array[x.col * i + j] = x.array[i][j] - y.array[i][j];
			}
		}

		return Matrix(temp_array, x.row, x.col);
	}

	else {

		std::cout << "型が正しくありません." << std::endl;
		return x.Zero(1, 1);
	}
}

//pro関数
template<class Type>
Matrix<Type> Matrix<Type>::pro(const Matrix& x, const Matrix& y) const {

	Type* temp_array = nullptr;
	temp_array = new Type[y.row * x.col];

	if (x.col == y.row) {

		for (int i = 0; i < x.col; i++)
			for (int j = 0; j < y.row; j++) {

				temp_array[i * x.col + j] = Sigma_array(x, y, i, j);
			}

		return Matrix(temp_array, y.row, x.col);
	}

	else {

		std::cout << "型が正しくありません." << std::endl;
		return x.Zero(1, 1);
	}
}

//pro_n関数
template<class Type>
Matrix<Type> Matrix<Type>::pro_n(const Matrix& x, Type nd) const {

	Type* temp_array = nullptr;
	temp_array = new Type[x.row * x.col];

	for (int i = 0; i < x.row; i++) {
		for (int j = 0; j < x.col; j++) {

			temp_array[i * x.col + j] = (x.array[i][j]) * nd;
		}
	}

	return Matrix(temp_array, x.row, x.col);
}


//ゼロ行列を生成
template<class Type>
Matrix<Type> Matrix<Type>::Zero(int _row, int _col) const {

	Type* temp_array = nullptr;
	temp_array = new Type[_row * _col];

	for (int i = 0; i < _row; i++) {
		for (int j = 0; j < _col; j++) {

			temp_array[_col * i + j] = 0;
		}
	}

	return Matrix(temp_array, _row, _col);
}


//単位行列を生成
template<class Type>
Matrix<Type> Matrix<Type>::Identity(int _row, int _col) const {

	if (_row == _col) {

		Type* temp_array = nullptr;
		temp_array = new Type[_row * _col];

		for (int i = 0; i < _row; i++)
			for (int j = 0; j < _col; j++) {

				if (i == j)
					temp_array[_col * i + j] = 1;
				else
					temp_array[_col * i + j] = 0;
			}

		return Matrix(temp_array, _row, _col);
	}

	else {

		std::cout << "型が正しくありません." << std::endl;
		return Zero(1, 1);
	}
}


//転置行列にする
template<class Type> 
Matrix<Type> Matrix<Type>::T(){

	Type* temp_array = nullptr;
	temp_array = new Type[this->row * this->col];

	for (int i = 0; i < this->row; i++) {
		for (int j = 0; j < this->col; j++) {
			temp_array[this->col * j + i] = this->array[i][j];
		}
	}

	return Matrix(temp_array, this->col, this->row); //行と列を逆にわたす
}


//正則判定
template<class Type>
bool Matrix<Type>::Is_regular() {

	Type* temp_array = nullptr;
	temp_array = new Type[this->row * this->col];

	for (int i = 0; i < this->row; i++) {
		for (int j = 0; j < this->col; j++) {

			temp_array[i * this->col + j] = this->array[i][j];
		}
	}

	Determinant<Type> A(temp_array, this->row, this->col);
	return (bool)A.det();
}


//余因子を求める
template<class Type>
Type Matrix<Type>::cof(int _i, int _j) {
	//_i行, _j列の余因子を返す
	//デフォルトは1行１列で展開

	if (this->get_row() != this->get_col()) {

		std::cout << "この関数は、正方行列のみに対応しています." << std::endl;

		return -1;
	}

	int num = 0;

	Type* temp_array = nullptr;
	temp_array = new Type[(this->row - 1) * (this->col - 1)];

	int bit = ((_i + 1) + (_j + 1)) % 2 ? -1 : 1;

	for (int i = 0; i < this->row; i++) {

		if (i == _i) continue;
		for (int j = 0; j < this->col; j++) {

			if (j == _j) continue;

			temp_array[num] = this->array[i][j];
			num++;
		}
	}

	Determinant<Type> A(temp_array, this->row - 1, this->col - 1);
	return A.det() * bit;
}


//逆行列を返す
template<class Type> 
Matrix<Type> Matrix<Type>::Inverse() {

	if (Is_regular()) {
		//正則ならば
		Type* temp_array = nullptr;
		temp_array = new Type[this->row * this->col];

		for (int i = 0; i < this->row; i++) {
			for (int j = 0; j < this->col; j++) {

				temp_array[i * this->col + j] = this->array[i][j];
			}
		}

		Determinant<Type> A(temp_array, this->row, this->col);

		return this->cof_matrix() * (1 / A.det());
	}
	else {

		std::cout << "この関数は正則な行列にのみ対応しています." << std::endl;
		return *this;
	}
}


//余因子行列を返す
template<class Type>
Matrix<Type> Matrix<Type>::cof_matrix() {

	Matrix<Type> A = A.Zero(this->row, this->col);
	for (int i = 0; i < this->row; i++) {
		for (int j = 0; j < this->col; j++) {

			A.array[i][j] = cof(i, j);
		}
	}

	Matrix<Type> B = A.T();
	std::cout << B << std::endl;
	return A.T(); //転置してます
}


//階級を求める
template<class Type> 
int Matrix<Type>::Rank() {

	int rank = 0;

	this->Sweep(); //掃き出しを行う

	for (int i = 0; i < this->row; i++) {
		
		for (int j = 0; j < this->col; j++) {

			if (this->array[i][j] != 0) {
				
				rank++;
				break;
			}
		}
	}

	return rank;
}


