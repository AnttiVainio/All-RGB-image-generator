#include <fstream>

inline void updateHeader(char *headerPos, const unsigned int value) {
	headerPos[0] = value & 0xff;
	headerPos[1] = (value >> 8) & 0xff;
	headerPos[2] = (value >> 16) & 0xff;
	headerPos[3] = value >> 24;
}

bool saveBMP(const unsigned char *data, const char *filepath, const unsigned int width, const unsigned int height) {
	std::ofstream file(filepath, std::ios::binary);
	if (!file.good()) {
		file.close();
		return false;
	}

	// Create the header
	char header[54] = {
		66, 77,      // BM
		0, 0, 0, 0,  // size of the file
		0, 0, 0, 0,
		54, 0, 0, 0, // offset to image data
		40, 0, 0, 0, // size of this header
		0, 0, 0, 0,  // width of the image
		0, 0, 0, 0,  // height of the image
		1, 0, 24, 0, // bits per pixel
		0, 0, 0, 0,
		0, 0, 0, 0,  // size of the pixel data
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0
	};

	const unsigned char padding = width % 4;
	const char *paddingArray = header + 50; // use an empty part of the header for writing padding to the file

	// Update the header
	updateHeader(header + 2, width * height * 3 + padding * height + 54); // size of the file
	updateHeader(header + 18, width); // width of the image
	updateHeader(header + 22, height); // height of the image
	updateHeader(header + 34, width * height * 3); // size of the pixel data
	file.write(header, 54);

	// Save the pixel data
	const unsigned char *dataPos = data;
	for(unsigned int i = 0; i < height; i++) {
		for(unsigned int j = 0; j < width; j++) {
			file.put(dataPos[2]);
			file.put(dataPos[1]);
			file.put(dataPos[0]);
			dataPos += 3;
		}
		file.write(paddingArray, padding);
	}

	file.close();
	return true;
}
