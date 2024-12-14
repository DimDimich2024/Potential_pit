
// Potential_pit.h: главный файл заголовка для приложения PROJECT_NAME
//

#pragma once

#ifndef __AFXWIN_H__
	#error "включить pch.h до включения этого файла в PCH"
#endif

#include "resource.h"		// основные символы


// CPotentialpitApp:
// Сведения о реализации этого класса: Potential_pit.cpp
//

class CPotentialpitApp : public CWinApp
{
private:
	Gdiplus::GdiplusStartupInput input;
	ULONG_PTR token;
public:
	CPotentialpitApp();

// Переопределение
public:
	virtual BOOL InitInstance();
	virtual BOOL ExitInstance();

// Реализация

	DECLARE_MESSAGE_MAP()
};

extern CPotentialpitApp theApp;
