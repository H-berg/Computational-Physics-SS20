#include <fstream>
#include <vector>
#include <Eigen/Dense>


int loadData(Eigen::MatrixXd &mat, std::string filename, size_t row, size_t col)
{

	std::ifstream input;
	input.open(filename, std::fstream::in | std::fstream::binary);
	if (input.is_open())
	{
		std::vector<unsigned char> vec((std::istreambuf_iterator<char>(input)),
				std::istreambuf_iterator<char>() );
                mat= Eigen::Map<Eigen::Matrix<unsigned char, Eigen::Dynamic,Eigen::Dynamic> >(vec.data(),row,col).template cast<double>();
		input.close();
		return 1;
	}
	return 0;

}
int storeData(Eigen::MatrixXd &mat, std::string filename)
{
	Eigen::Matrix<char, Eigen::Dynamic, Eigen::Dynamic> pic= mat.template cast<char>();
	std::ofstream output;
	output.open(filename,std::fstream::trunc);
	if (output.is_open())
	{
		output.write(pic.data(),pic.size());
		output.close();
		return 1;
	}
	return 0;


}
