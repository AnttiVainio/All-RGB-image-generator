#include <iostream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#if defined(_OPENMP)
	#include <omp.h>
#endif

// This might not be defined in cmath
#define M_PI 3.14159265358979323846

/*
 * NOTE:
 * There are also some OpenMP pragmas in the code that will speed up the program when compiled with OpenMP
 * However, if enabled, the program will not go through the coordinates in a completely deterministic order anymore
 *   but otherwise the program works correctly and the images will contain only unique colors
 */

// Using this SIZE will make sure that the program uses exactly all colors in the generated color pool
// Setting arbitrary IMG_WIDTH and IMG_HEIGHT may give the program more colors to work with and make it easier to create smooth images
#define SIZE 8 // final image size in pixels will be this number cubed (max 16)
#define IMG_WIDTH (SIZE * SIZE * SIZE)
#define IMG_HEIGHT (SIZE * SIZE * SIZE)

#define EXTRA_COLORS 0 // can be used to increase the unique color pool to produce smoother images
#define START_AMOUNT 100 // amount of initial random pixels
#define FIND_AMOUNT 20 // amount of pixels to find before choosing color
#define FIND_SIZE 100 // distance to search for pixels before aborting
#define COLOR_SEARCH_RANGE 100 // amount of colors to check for a new pixel from the remaining color pool
#define COLOR_ACCEPTANCE 25 // absolute acceptable difference between r, g and b values

#define SAVE_INTERVAL 0.1 // progress between image saves (max progress = 1.0)

// This type setting affects the shape of the generated image by generating the pixels in different order
// 0: top to bottom
// 1: random
// 2: rotate around the middle
// 3: middle to edges
// 4: edges to middle
// 5: multiple circles
// 6: another multiple circles
// 7: weirdness
// 8: vertical lines
// 9: random spots
#define TYPE 6

// Setting some random range here will affect the sorting of the different pixel positions
// Good values depend on the selected TYPE and usually vary from (-1 to 1) to (-0.001 to 0.001)
#define RANDOM_ORDER_MIN -0.1
#define RANDOM_ORDER_MAX 0.1

// A class for wrapping some C++ random functions
class RNG {
	private:
		std::default_random_engine rng;
	public:
		RNG(): rng(std::chrono::system_clock::now().time_since_epoch().count()) {}
		template <typename T> T getRandom(const T a, const T b) {
			return std::uniform_int_distribution<T>(a, b)(rng);
		}
		template <typename T> T getRandomf(const T a, const T b) {
			return std::uniform_real_distribution<T>(a, b)(rng);
		}
		std::default_random_engine::result_type operator()() {
			return rng();
		}
};

bool saveBMP(const unsigned char *data, const char *filepath, const unsigned int width, const unsigned int height);

typedef std::pair<unsigned int, float> Color;
typedef std::pair<unsigned int, float> Coordinate;

bool compareSecond(const std::pair<unsigned int, float> &a, const std::pair<unsigned int, float> &b) {
	return a.second < b.second;
}

float getHue(const unsigned int c) {
	return atan2(sqrt(3.0) * (((c >> 8) & 255) - ((c >> 16) & 255)), 2.0 * (c & 255) - ((c >> 8) & 255) - ((c >> 16) & 255));
}

int main() {
	const auto benchmark = std::chrono::high_resolution_clock::now();

	#if defined(_OPENMP)
		std::cout << "Using " << omp_get_max_threads() << " threads" << std::endl << std::endl;
	#endif

	RNG rng;

	unsigned int COLOR_SIZE = 1; // final amount of unique colors will be this number cubed
	while (COLOR_SIZE * COLOR_SIZE * COLOR_SIZE < IMG_HEIGHT * IMG_WIDTH) COLOR_SIZE += 1;
	COLOR_SIZE = std::min(COLOR_SIZE + EXTRA_COLORS, 256u); // add some extra colors for smoother images

	const int imgNameLen = (int)log10(1.0 / SAVE_INTERVAL) + 1;

	const std::unique_ptr<unsigned char[]> data(new unsigned char[IMG_WIDTH * IMG_HEIGHT * 3]());

	std::cout << "Image size: " << IMG_WIDTH << " " << IMG_HEIGHT << std::endl;
	std::cout << "Color channel size: " << COLOR_SIZE << std::endl;
	std::cout << (COLOR_SIZE * COLOR_SIZE * COLOR_SIZE - IMG_HEIGHT * IMG_WIDTH) << " extra colors" << std::endl << std::endl;

	std::cout << "Generating colors and coordinates" << std::endl;

	std::vector<Coordinate> coords(IMG_WIDTH * IMG_HEIGHT);
	std::vector<Color> colors;

	// Generate random spots for TYPE 9
	std::vector<std::pair<float, float>> randomSpots(50);
	for (auto i = randomSpots.begin(); i != randomSpots.end(); i++) {
		i->first = rng.getRandomf<float>(0, 1);
		i->second = rng.getRandomf<float>(0, 1);
	}

	#pragma omp parallel
	{
		#pragma omp single nowait
		{
			// Generate color values for a single channel
			std::vector<unsigned int> cvals(COLOR_SIZE);
			for (unsigned int i = 0; i < COLOR_SIZE; i++) {
				cvals[i] = (unsigned int)(round(i * 255.0 / float(COLOR_SIZE - 1)));
			}

			// Generate colors
			for (auto r = cvals.begin(); r != cvals.end(); r++) {
				for (auto g = cvals.begin(); g != cvals.end(); g++) {
					for (auto b = cvals.begin(); b != cvals.end(); b++) {
						const unsigned int color = *r + (*g << 8) + (*b << 16);
						colors.push_back(Color(color, getHue(color)));
					}
				}
			}
		}

		// Generate coordinates
		#pragma omp for
		for (unsigned int i = 0; i < IMG_WIDTH * IMG_HEIGHT; i++) {
			coords[i].first = i;
			// Create some sorting values for coordinates
			const float x = float(i % IMG_WIDTH) / IMG_WIDTH;
			const float y = float(i / IMG_WIDTH) / IMG_HEIGHT;
			const float dist = sqrt((0.5 - x) * (0.5 - x) + (0.5 - y) * (0.5 - y)); // distance to middle

			switch (TYPE) {
				case 0: // top to bottom
					coords[i].second = -i;
					break;
				case 1: // random
					coords[i].second = rng();
					break;
				case 2: // rotate around the middle
					coords[i].second = std::atan2(x - 0.5, y - 0.5);
					break;
				case 3: // middle to edges
					coords[i].second = dist;
					break;
				case 4: // edges to middle
					coords[i].second = -dist;
					break;
				case 5: // multiple circles
					coords[i].second = sin(((0.5 - x) * (0.5 - x) + (0.5 - y) * (0.5 - y)) * 12.0 * M_PI);
					break;
				case 6: // another multiple circles
					coords[i].second = sin(100.0 * -dist) + dist * 15.0 + fmod(atan2(x - 0.5, y - 0.5) + dist * 10.0 + M_PI + rng.getRandomf<float>(-0.1, 0.1), (2.0 * M_PI));
					break;
				case 7: // weirdness
					coords[i].second = std::pow(i / IMG_WIDTH - (cos(i * 39.999 / IMG_HEIGHT) + 0.5) * IMG_HEIGHT, 4)
					                 + std::pow(i % IMG_WIDTH - (sin(i * 40.0 / IMG_WIDTH) + 0.5) * IMG_WIDTH, 4);
					break;
				case 8: // vertical lines
					coords[i].second = i + (i % IMG_WIDTH % 10) * IMG_WIDTH * IMG_HEIGHT;
					break;
				case 9: // random spots
					coords[i].second = -2;
					for (auto j = randomSpots.cbegin(); j != randomSpots.cend(); j++) {
						coords[i].second = std::max(coords[i].second, -(j->first - x) * (j->first - x) - (j->second - y) * (j->second - y));
					}
					break;
			}

			// Add some randomness
			if (RANDOM_ORDER_MIN && RANDOM_ORDER_MAX) {
				coords[i].second += rng.getRandomf<float>(RANDOM_ORDER_MIN, RANDOM_ORDER_MAX);
			}
		}
	}

	// These iterators mark the start of the remaining coordinates and colors
	auto coordIt = coords.begin();
	auto colorIt = colors.begin();

	// Create some random initial pixels
	std::cout << "Creating some random pixels" << std::endl;
	for (unsigned int i = 0; i < START_AMOUNT && i < IMG_WIDTH * IMG_HEIGHT; i++) {
		unsigned int r = rng.getRandom<unsigned int>(0, std::distance(coordIt, coords.end()) - 1);
		std::swap(*coordIt, *(coordIt + r));
		const unsigned int coord = coordIt->first;
		coordIt++;
		r = rng.getRandom<unsigned int>(0, std::distance(colorIt, colors.end()) - 1);
		std::swap(*colorIt, *(colorIt + r));
		const unsigned int color = colorIt->first;
		colorIt++;
		data[coord * 3    ] = color & 255;
		data[coord * 3 + 1] = (color >> 8) & 255;
		data[coord * 3 + 2] = (color >> 16) & 255;
	}

	std::cout << "Sorting colors and coordinates" << std::endl;

	#pragma omp parallel
	{
		#pragma omp single
		{
			// Sort the colors by their hue
			#pragma omp task
			std::sort(colorIt, colors.end(), compareSecond);
			// Sort coordinates based on the sorting index
			std::sort(coordIt, coords.end(), compareSecond);
		}
	}

	std::cout << "Creating image" << std::endl << std::endl;

	auto starttime = std::chrono::high_resolution_clock::now();
	unsigned int total = 0;
	unsigned int randomed = 0;
	unsigned int acceptable = 0;
	float save = SAVE_INTERVAL;
	unsigned int imgName = 0;

	// Generate the image now
	const unsigned int cend = coords.size();
	// Dynamic schedule is actually quite fast and generates images more consistently (coordinates closer to the correct order)
	#pragma omp parallel for schedule(dynamic)
	for (unsigned int c = START_AMOUNT; c < cend; c++) {
		std::ostringstream filename;
		bool saveImage = false;

		#pragma omp critical(progress)
		{
			const double progress = total * 100.0 / IMG_WIDTH / IMG_HEIGHT;
			total++;

			// Check whether to save the image
			if (progress / 100.0 > save) {
				filename << "images/" << std::setfill('0') << std::setw(imgNameLen) << imgName << ".bmp";
				save += SAVE_INTERVAL;
				imgName++;
				saveImage = true;
			}

			// Print progress
			const auto endtime = std::chrono::high_resolution_clock::now();
			if (c == coords.size() - 1 || std::chrono::duration_cast<std::chrono::seconds>(endtime-starttime).count() > 0) {
				starttime = endtime;
				std::cout << std::fixed << std::setprecision(1)
					<< progress << "%  "
					<< (randomed * 100.0 / total) << "% randomed  "
					<< (acceptable * 100.0 / (total - randomed)) << "% acceptable" << std::endl;
			}
		}

		// Save bmp
		if (saveImage) saveBMP(data.get(), filename.str().c_str(), IMG_WIDTH, IMG_HEIGHT);

		const unsigned int coord = coords[c].first;
		const int x = coord % IMG_WIDTH;
		const int y = coord / IMG_WIDTH;
		unsigned int found = 0;
		float r = 0, g = 0, b = 0;
		int i = 1;
		unsigned int color;

		// Search other pixels near the current coordinate
		while (found < FIND_AMOUNT) {
			// Go through the different directions in a random order
			unsigned char randOrder[] = {0,1,2,3};
			std::random_shuffle(randOrder, randOrder + 4);
			for (unsigned char *k = randOrder; k < randOrder + 4; k++) {
				if ((*k == 0 && x - i >= 0) || (*k == 1 && x + i < IMG_WIDTH)) {
					for (int j = -i; found < FIND_AMOUNT && j <= i; j++) {
						const unsigned int pos = ((y + j) * IMG_WIDTH + x + (*k == 0 ? -i : i)) * 3;
						if (y + j >= 0 && y + j < IMG_HEIGHT && data[pos] + data[pos + 1] + data[pos + 2]) {
							found++;
							r += data[pos];
							g += data[pos + 1];
							b += data[pos + 2];
						}
					}
				}
				else if ((*k == 2 && y - i >= 0) || (*k == 3 && y + i < IMG_HEIGHT)) {
					for (int j = -i + 1; found < FIND_AMOUNT && j < i; j++) {
						const unsigned int pos = ((y + (*k == 2 ? -i : i)) * IMG_WIDTH + x + j) * 3;
						if (x + j >= 0 && x + j < IMG_WIDTH && data[pos] + data[pos + 1] + data[pos + 2]) {
							found++;
							r += data[pos];
							g += data[pos + 1];
							b += data[pos + 2];
						}
					}
				}
			}
			i++;
			if (i == FIND_SIZE) found = FIND_AMOUNT + 1;
		}

		// Not enough pixels were found before aborting - choosing random color
		if (found == FIND_AMOUNT + 1) {
			#pragma omp critical(colorSelection)
			{
				const unsigned int r = rng.getRandom<unsigned int>(0, std::distance(colorIt, colors.end()) - 1);
				std::swap(*colorIt, *(colorIt + r)); // swapping colors like this makes the results more noisy but speeds up the program
				color = colorIt->first;
				colorIt++;
				randomed++;
			}
		}
		// Enough pixel were found
		else {
			r /= FIND_AMOUNT;
			g /= FIND_AMOUNT;
			b /= FIND_AMOUNT;
			// Look at the color vector based on the hue of the desired color
			const float tempPos = (atan2(sqrt(3.0) * (g - b), 2.0 * r - g - b) + M_PI) / M_PI / 2.0001;
			int i = 0;
			float best = 1000;
			unsigned int stop = 0;
			#pragma omp critical(colorSelection)
			{
				const int clen = std::distance(colorIt, colors.end());
				const int pos = int(tempPos * clen);
				int index = pos;
				// Look for a good color
				while (stop < COLOR_SEARCH_RANGE && best > COLOR_ACCEPTANCE) {
					const int i2 = ((i + pos) % clen + clen) % clen;
					const unsigned int color2 = (colorIt + i2)->first;
					const float dist = std::max(std::max(
						std::fabs(r - (color2 & 255)),
						std::fabs(g - ((color2 >> 8) & 255))),
						std::fabs(b - ((color2 >> 16) & 255)));
					if (dist < best) {
						best = dist;
						index = i2;
					}
					i = i * -2 + 1; // get the next index to be looked at around the previous indexes
					stop++;
				}
				if (best <= COLOR_ACCEPTANCE) acceptable++;
				std::swap(*colorIt, *(colorIt + index)); // swapping colors like this makes the results more noisy but speeds up the program
				color = colorIt->first;
				colorIt++;
			}
		}

		data[coord * 3    ] = color & 255;
		data[coord * 3 + 1] = (color >> 8) & 255;
		data[coord * 3 + 2] = (color >> 16) & 255;
	}

	std::ostringstream filename;
	filename << "images/" << std::setfill('0') << std::setw(imgNameLen) << imgName << ".bmp";
	saveBMP(data.get(), filename.str().c_str(), IMG_WIDTH, IMG_HEIGHT);

	const auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::endl << "Time taken: "
		<< ((std::chrono::duration_cast<std::chrono::nanoseconds>(end-benchmark).count() / 10000000) / 100.0)
		<< " seconds" << std::endl;
	return 0;
}
