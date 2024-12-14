
// Potential_pitDlg.h: файл заголовка
//

#pragma once
#include <vector>
#include <cmath>


// Диалоговое окно CPotentialpitDlg
class CPotentialpitDlg : public CDialogEx
{
// Создание
public:
	CPotentialpitDlg(CWnd* pParent = nullptr);	// стандартный конструктор

// Данные диалогового окна
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_POTENTIAL_PIT_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// поддержка DDX/DDV


// Реализация
protected:
	HICON m_hIcon;

	// Созданные функции схемы сообщений
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()


public:
	int n = 1;
	int k; //уровень энергии
	int k_max=10; //максимальное количество уровней энергии
	double ek; //значение энергии на уровне k
	double R;
	double k_potential; //наклон прямой
	std::vector <long double> phi; //фаза фи
	std::vector <long double> energy; //энергия
	std::vector <long double> wave; //волновая функция
	std::vector <long double> z; //координата

	//фазовая функция
	void DrawPhi(std::vector <long double> phi, std::vector <long double> e);
	void DrawWave(std::vector <long double> wave, std::vector <long double> z);
	afx_msg void OnBnClickedOk();
	
	//double C;
	double C;
};



