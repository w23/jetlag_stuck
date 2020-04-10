#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <windows.h>

#define COUNTOF(a) (sizeof((a))/sizeof((a)[0]))

//#define DEFAULT_FONT L"Tahoma"
#define DEFAULT_FONT L"Consolas"
// L"Arial"

struct Label {
	//const char *uniform;
	const char *ansi;
	//const wchar_t *unicode;
	//const wchar_t *font;
};

//const wchar_t rojikoma[] = { 0x30ed, 0x30b8, 0x30b3, 0x30de, 0x0020, 0x0425, 0x0423, 0x0419, 0 };

Label labels[] = {
	"Notice me senpai\n"
	"Farbrausch\n"
	"Logicoma\n"
	"LJ\nPrismbeings\nAlcatraz\nConspiracy\nQuite\nMercury\nSandS\nTitan\nThrob\n"
};

const unsigned int TEXT_WIDTH = 1024, TEXT_HEIGHT = 1024;
static BITMAPINFO bitmap_info = {{
	/*biSize = */sizeof(BITMAPINFOHEADER),
	/*biWidth = */TEXT_WIDTH,
	/*biHeight = */TEXT_HEIGHT,
	/*biPlanes = */1,
	/*biBitCount = */32,
	/*biCompression = */BI_RGB,
	0,
	0,
	0,
	0,
	0,
}, {0}};

int main() {
	void *bitmap_ptr = NULL;
	const HDC text_dc = CreateCompatibleDC(NULL);
	const HBITMAP dib = CreateDIBSection(text_dc, &bitmap_info, DIB_RGB_COLORS, &bitmap_ptr, NULL, 0);
	const HGDIOBJ obj = SelectObject(text_dc, dib);
	RECT rect = { 0, 0, TEXT_WIDTH, TEXT_HEIGHT };
	SetTextColor(text_dc, RGB(255, 255, 255));
	SetBkMode(text_dc, TRANSPARENT);

	int size = 44;
	//HFONT font = CreateFont(size, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*CLEARTYPE_QUALITY*/ NONANTIALIASED_QUALITY, 0, DEFAULT_FONT);
	HFONT font = CreateFont(size, 0, 0, 0, 700, 0, 0, 0, 0, 0, 0, /*CLEARTYPE_QUALITY*/ ANTIALIASED_QUALITY, 0, DEFAULT_FONT);
	SelectObject(text_dc, font);
	for (int i = 0; i < COUNTOF(labels); ++i) {
		const Label &l = labels[i];

		// needs rect; keeps newlines
		if (l.ansi) {
			DrawTextA(text_dc, l.ansi, -1, &rect, DT_CALCRECT);
			DrawTextA(text_dc, l.ansi, -1, &rect, 0);
			//DrawTextA(text_dc, l.ansi, -1, &rect, DT_CENTER | DT_VCENTER | DT_CALCRECT);
			//DrawTextA(text_dc, l.ansi, -1, &rect, DT_CENTER | DT_VCENTER);
			// DT_SINGLELINE | DT_NOCLIP);
		/*} else {
			DrawTextW(text_dc, l.unicode, -1, &rect, DT_CALCRECT);
			DrawTextW(text_dc, l.unicode, -1, &rect, 0);
		*/
		}

		/*
		l.x = rect.left;
		l.y = rect.top; l.w = rect.right - rect.left;
		l.h = rect.bottom - rect.top;
		*/
		printf("%s %d %d %d %d\n", l.ansi, rect.left, rect.top, rect.right, rect.bottom);

		rect.left = 0;
		rect.right = TEXT_WIDTH;
		rect.top = rect.bottom;
		rect.bottom = TEXT_HEIGHT;

		/*
		// ignores newlines
		TextOutW(text_dc, 0, 400, ro, 8);
		TextOutA(text_dc, 0, 600, shader_program, strlen(shader_program));
		*/
	}

	//stbi_flip_vertically_on_write(1);
	int result = stbi_write_png("texture.png", TEXT_WIDTH, TEXT_HEIGHT, 4, bitmap_ptr, 4 * TEXT_WIDTH);
	printf("stbi result = %d\n", result);
	result = stbi_write_jpg("texture.jpg", TEXT_WIDTH, TEXT_HEIGHT, 4, bitmap_ptr, 80);
	printf("stbi result = %d\n", result);

	DeleteObject(font);
	DeleteObject(dib);
	DeleteDC(text_dc);
	return 0;
}
