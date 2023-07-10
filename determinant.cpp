#include"Matrix.h"


template class Determinant<double>;
template class Determinant<int>;

template<class Type>
void Determinant<Type>::set_point(int point2row, int point2col) {
	//着目対象の行または列

	if ((point2row == -1 || point2col == -1) && (point2row >= 0 || point2col >= 0)) {

		pointed_row = point2row;
		pointed_col = point2col;
	}

	else {

		std::cout << "指定方式エラー" << std::endl;
	}
}

//要素がすべてゼロの行列式を生成
template<class Type>
Determinant<Type> Determinant<Type>::Zero(int _row, int _col) const {

	Type* temp_array = nullptr;
	temp_array = new Type[_row * _col];

	for (int i = 0; i < _row; i++) {
		for (int j = 0; j < _col; j++) {

			temp_array[_col * i + j] = 0;
		}
	}

	return Determinant(temp_array, _row, _col);
}

//要素が単位行列の行列式を生成
template<class Type>
Determinant<Type> Determinant<Type>::Identity(int _row, int _col) const {

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

		return Determinant(temp_array, _row, _col);
	}

	else {

		std::cout << "型が正しくありません." << std::endl;
		return Zero(1, 1);
	}
}


//行列式を求める
template<class Type>
Type Determinant<Type>::det() const {

	if (this->get_row() != this->get_col()) {

		std::cout << "この関数は、正方行列のみに対応しています." << std::endl;
		return -1;
	}

	Type V = 0;


	if (this->get_row() == 1) {

		return this->array[0][0];
	}

	else if (this->get_row() > 2) {

		for (int n = 0; n < this->row; n++) {

			Determinant A = cof_for_det(0, n);
			V += A.det();
		}
	}

	else if (this->get_row() == 2) {

		V += this->array[0][0] * this->array[1][1] - this->array[0][1] * this->array[1][0];
	}


	return V;
}


//余因子を求める(行列式計算用)_privateメンバ
template<class Type>
Determinant<Type> Determinant<Type>::cof_for_det(int _i, int _j) const {

	//デフォルトは1行１列で展開


	int num = 0;

	Type P = this->array[_i][_j];

	int bit = ((_i + 1) + (_j + 1)) % 2 ? -1 : 1;

	Type* temp_array = nullptr;
	temp_array = new Type[(this->row - 1) * (this->col - 1)]; //n-1次に下げた箱

	for (int i = 0; i < this->row; i++) {

		if (i == _i) continue;
		for (int j = 0; j < this->col; j++) {

			if (j == _j) continue;


			//行列式の乗算を考慮
			if (_i != 0) {
				if (i == 0) {

					temp_array[num] = this->array[i][j] * bit * P;
					num++;
				}
				else {

					temp_array[num] = this->array[i][j];
					num++;
				}
			}

			else {

				if (i == 1) {

					temp_array[num] = this->array[i][j] * bit * P;
					num++;
				}
				else {

					temp_array[num] = this->array[i][j];
					num++;
				}
			}

		}
	}

	return Determinant(temp_array, this->row - 1, this->col - 1);
}
