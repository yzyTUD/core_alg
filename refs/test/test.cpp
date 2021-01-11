
#include <zrlib.h>	//

int main(){

    /// testing logic lib.
		(new zrGraph())->gen_test_data_and_test();
		(new zrTree())->gen_test_data_and_test();
		zr_hannoi(2, data_zr_hannoi_1[0], data_zr_hannoi_1[1], data_zr_hannoi_1[2]);

	/// testing math lib.
		zr_calcu_num_of_rabbits(5);
		zr_calcu_mandelbrot_output_as_ppm();
		zr_calcu_series();

	/// testing util lib.
		std::cout << zr_get_current_time() << std::endl;
		zr_double_the_given_file("1.txt");
		zr_bitwise_reverse("1.txt");
}