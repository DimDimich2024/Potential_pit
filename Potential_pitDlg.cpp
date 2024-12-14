
// Potential_pitDlg.cpp: файл реализации
//

#include "pch.h"
#include "framework.h"
#include <algorithm>
#include <fstream>
#include <thread>
#include "Potential_pit.h"
#include "Potential_pitDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// Диалоговое окно CAboutDlg используется для описания сведений о приложении
double random(double min, double max) {
	return min + (max - min) * rand() / RAND_MAX;
}
const long double pi = 3.14159265358979323846;

long double constV0 = 8., constA = 5.;
int NN = 0, last_k_max = 0;
long double UU = 0, wwidth_between = 0;
template <typename T>
int sign(T x) {
	if (x > 0) return 1;
	else if (x == 0) return 0;
	else return -1;
}

long double Pow2(long double x)
{
	return x * x;
}
long double PotentialFunction(long double z, long double k, double C)
{
	//ширина ямы = доля от всей ширины / кол-во ям
	if (z <= 0)
		return C;
	else
		return z * k;
}

long double PhaseFunction(long double fi, long double z, long double e, long double k, double C)
{
	long double potential = PotentialFunction(z, k,C);
	return (Pow2(cos(fi)) * (potential - e) - Pow2(sin(fi)));
}
long double RadiusFunction(long double r, long double e, long double z, long double fi, long double k, double C)
{
	long double potential = PotentialFunction(z, k,C);
	return r * (1 + potential - e) * cos(fi) * sin(fi);
}
long double SolutionPF(int n)
{
	return -0.5 * pi * (2. * n + 1.);
}

//возвращает фи в R
long double MethodRK4(long double R, long double e, long double k, double C)
{
	int numStep = 500;
	long double step = R / numStep;

	std::vector<long double> massR(numStep + 1);
	std::vector<long double> massFi(numStep + 1);
	massFi[0] = pi / 2.;
	massR[0] = 0;

	for (int i = 0; i < numStep; i++)
	{
		long double r = step*i;

		long double k1 = PhaseFunction(massFi[i], r, e, k,C);
		long double k2 = PhaseFunction(massFi[i] + step * k1 / 2, r + step / 2, e, k,C);
		long double k3 = PhaseFunction(massFi[i] + step * k2 / 2, r + step / 2, e, k,C);
		long double k4 = PhaseFunction(massFi[i] + step * k3, r + step, e, k,C);

		massFi[i + 1] = massFi[i] + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		massR[i + 1] = step * (i+1);
	}

	return massFi[numStep];
}

// фи от е в след точке
long double MethodRK4(long double lastfi, long double z, long double e, long double step, long double k, double C)
{

	long double k1 = PhaseFunction(lastfi, z, e, k,C);
	long double k2 = PhaseFunction(lastfi + step * k1 / 2, z + step / 2, e, k,C);
	long double k3 = PhaseFunction(lastfi + step * k2 / 2, z + step / 2, e, k,C);
	long double k4 = PhaseFunction(lastfi + step * k3, z + step, e, k,C);

	return lastfi + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

std::vector<long double> CreateWaveFunction(long double R, long double ek, std::vector<long double>& masZ, long double fi_last, long double k, double C) //masZ меняется
{
	int numStep = 5000;
	long double step = R / numStep;

	std::vector<long double> masR(numStep + 1);
	std::vector<long double> masPsi(numStep + 1);
	masZ.clear();
	masZ = std::vector<long double>(numStep + 1,0);
	//std::vector<long double>phi;
	//MethodRK4(R, ek, phi);
	masR[0] =1;
	masZ[0] = -R;
	masPsi[0] = 0;
	long double fik = pi / 2;
	std::fstream fout("WaveFunc.txt", std::ios::out);
	fout << "r\tphi\tz" << std::endl;
	for (int i = 0; i < numStep; i++)
	{
		long double zi =-R+2* i*step;
	
		//long double fik = MethodRK4(zi + step, ek);
		double fik1 = MethodRK4(fik, zi, ek, step, k,C);
		
		
		
		fout << masPsi[i] << '\t'<< fik << '\t' << zi << std::endl;
		long double k1 = RadiusFunction(masR[i], ek, zi, fik, k,C);
		long double k2 = RadiusFunction(masR[i] + step * k1 / 2, ek, zi + step / 2, (fik + fik1) / 2, k,C);
		long double k3 = RadiusFunction(masR[i] + step * k2 / 2, ek, zi + step / 2, (fik + fik1) / 2, k,C);
		long double k4 = RadiusFunction(masR[i] + step * k3, ek, zi + step, fik1, k,C);
		masR[i + 1] = masR[i] + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
		masPsi[i + 1] = masR[i + 1] * cos(fik1);
		masZ[i + 1] =-R+ 2*step*(i+1);
		fik = fik1;
	}

	fout.close();
	return masPsi;
}

std::vector<long double> CreatePhaseFunction(long double R, int k, std::vector<long double>& masE, long double k_potential, double C)//masE меняется
{
	long double solFi = SolutionPF(k + 1);

	long double step = 0.1/R; 	

	std::vector<long double> listFi;
	std::vector<long double> listE;

	for (int i = 0; ; i++)
	{
		long double e = step * i;
		long double valueFi = MethodRK4(R, e, k_potential,C);
		listFi.push_back(valueFi);
		listE.push_back(e);

		if (valueFi <= solFi) break;
	}

	masE = listE;
	return listFi;
}


long double MethodDevision(long double R, long double solFi, long double fi_right, long double e_right, long double& fik, long double k, double C) //fik - out
{
	//e_left, e_right - интервал оси х
	//
	long double accuracy = 1e-12;
	long double e_left = -1000, ek = 0;
	fik = 0;
	int iter = 0;
	double fi_left = MethodRK4(R, e_left, k,C);

	for (iter; iter < 10000; iter++)
	{
		//середина отрезка
		ek = (e_right + e_left) / 2;
		//значение фазы в этой точке
		fik = MethodRK4(R, ek, k,C);
		//если фаза в этой точке совпадает с нужной, прервем цикл
		if (abs(fik - solFi) <= accuracy)
		{
			break;
		}
		if ((fik - solFi)*(fi_right - solFi) > 0)
		{
			e_right = ek;
			//fi_right = fik;
			
		}
		else
		{
			e_left = ek;
		}
	}

	return ek;
}



class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// Данные диалогового окна
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // поддержка DDX/DDV

// Реализация
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// Диалоговое окно CPotentialpitDlg



CPotentialpitDlg::CPotentialpitDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_POTENTIAL_PIT_DIALOG, pParent)
	, k(0)
	, R(1)
	//, k_max(8)
	, ek(0)
	, k_potential(10)
	//, C(0)
	, C(15)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CPotentialpitDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_k, k);
	DDX_Text(pDX, IDC_R, R);
	DDX_Text(pDX, IDC_ek, ek);
	//	DDX_Text(pDX, IDC_k_max, k_max);
	DDX_Text(pDX, IDC_k_potential, k_potential);
	//DDX_Text(pDX, C, C);
	DDX_Text(pDX, IDC_C, C);
}

BEGIN_MESSAGE_MAP(CPotentialpitDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDOK, &CPotentialpitDlg::OnBnClickedOk)
END_MESSAGE_MAP()


// Обработчики сообщений CPotentialpitDlg

BOOL CPotentialpitDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Добавление пункта "О программе..." в системное меню.

	// IDM_ABOUTBOX должен быть в пределах системной команды.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != nullptr)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Задает значок для этого диалогового окна.  Среда делает это автоматически,
	//  если главное окно приложения не является диалоговым
	SetIcon(m_hIcon, TRUE);			// Крупный значок
	SetIcon(m_hIcon, FALSE);		// Мелкий значок

	// TODO: добавьте дополнительную инициализацию

	return TRUE;  // возврат значения TRUE, если фокус не передан элементу управления
}

void CPotentialpitDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// При добавлении кнопки свертывания в диалоговое окно нужно воспользоваться приведенным ниже кодом,
//  чтобы нарисовать значок.  Для приложений MFC, использующих модель документов или представлений,
//  это автоматически выполняется рабочей областью.

void CPotentialpitDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // контекст устройства для рисования

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Выравнивание значка по центру клиентского прямоугольника
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Нарисуйте значок
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		if (phi.size() > 0 || wave.size() > 0) {
			DrawPhi(phi, energy); //фи от энергии
			DrawWave(wave, z);
		}
		CDialogEx::OnPaint();
	}
}

// Система вызывает эту функцию для получения отображения курсора при перемещении
//  свернутого окна.
HCURSOR CPotentialpitDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CPotentialpitDlg::DrawPhi(std::vector <long double> phi, std::vector <long double> e)
{
	//рисование модели
	CDC* dc = GetDlgItem(IDC_phase)->GetDC();
	HDC hdc = *dc;
	Gdiplus::Graphics gr(hdc); //gr будет рисовать на IDC_Object
	CRect obj_rect; //прямоугольник области рисования
	GetDlgItem(IDC_phase)->GetClientRect(&obj_rect); //получаем размеры прямоугольника
	Gdiplus::Bitmap myBitmap(obj_rect.Width(), obj_rect.Height()); //создаем битовый образ
	
	Gdiplus::Graphics* grr = Gdiplus::Graphics::FromImage(&myBitmap); //создаем дополнительный объект класса для рисования объектов
	grr->SetSmoothingMode(Gdiplus::SmoothingModeHighSpeed); //устанавливаем сглаживание в режиме наилучшего качества
	
	grr->Clear(Gdiplus::Color::WhiteSmoke);//очистим фон
	

	Gdiplus::Matrix mtx(1.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f);//создаем матрицу с осью у, направленной вверх

	long double phi_max = *max_element(phi.begin(), phi.end()), 
			phi_min = *min_element(phi.begin(), phi.end()),
			e_max = *max_element(e.begin(), e.end());

	mtx.Translate(obj_rect.Width() * 0.05, obj_rect.Height() * (phi_max+ (0.04 * (phi_max - phi_min)))/ (phi_max - phi_min), Gdiplus::MatrixOrderAppend);//перенос начала координат
	mtx.Scale(obj_rect.Width() / (1.05*e_max), obj_rect.Height() / (1.04*(phi_max - phi_min))); //изменяем масштаб
	grr->SetTransform(&mtx); //применяем преобразования
	

	//рисуем сетку координат
	Gdiplus::Pen pGrid(Gdiplus::Color::Color(50, 50, 50), 0.01); //ручка для основной сетки
	Gdiplus::Pen pSubGrid(Gdiplus::Color::Color(200, 200, 200), 0.005); //ручка для основной сетки

	//рисуем оси
	grr->DrawLine(&pGrid, (Gdiplus::REAL)(-0.05 * e_max), (Gdiplus::REAL)(0), (Gdiplus::REAL)(e_max), (Gdiplus::REAL)(0));
	grr->DrawLine(&pGrid, (Gdiplus::REAL)(0), (Gdiplus::REAL)(phi_max + 0.04 * (phi_max - phi_min)), (Gdiplus::REAL)(0), (Gdiplus::REAL)(phi_min));

	//рисуем уровни энергий
	for (int i = 0; i <= k_max; i++) {
		int k = SolutionPF(i);
		grr->DrawLine(&pSubGrid, (Gdiplus::REAL)-0.1 * e_max, (Gdiplus::REAL)k, (Gdiplus::REAL)e_max, (Gdiplus::REAL)k);
	}

	Gdiplus::Pen pObj(Gdiplus::Color::Aqua, e_max / 500); //ручка для рисования объектов
	Gdiplus::SolidBrush brObj(Gdiplus::Color::Black); //ручка для рисования объектов
	//рисуем график
	int phi_size = phi.size();
	long double step = 0.0005;
	for (int i = 0; i < phi_size - phi_size*step;i+= phi_size*step) {
		grr->DrawLine(&pObj, (Gdiplus::REAL)e[i], (Gdiplus::REAL)phi[i], (Gdiplus::REAL)e[i+ phi_size*step], (Gdiplus::REAL)phi[i + phi_size*step]);
	}

	//сбрасываем матрицу преобразования
	mtx.Reset();
	mtx.Translate(obj_rect.Width() * 0.05, obj_rect.Height() * (phi_max + (0.04 * (phi_max - phi_min))) / (phi_max - phi_min), Gdiplus::MatrixOrderAppend);//перенос начала координат

	grr->SetTransform(&mtx); //применяем преобразования

	Gdiplus::REAL font_size = min((Gdiplus::REAL)(obj_rect.Width() / (2. * 20)), (Gdiplus::REAL)(obj_rect.Height() / (1. * 20)));
	//добавляем подписи
	Gdiplus::Font font(L"Arial", font_size, Gdiplus::FontStyleRegular, Gdiplus::UnitPixel); //создаем шрифт
	Gdiplus::SolidBrush brText(Gdiplus::Color::Black); //кисть для шрифта
	grr->SetTextRenderingHint(Gdiplus::TextRenderingHintClearTypeGridFit); //сглаживание в стиле cleartype

	//подпись оси x
	long double x = 0;
	for (long double i = -obj_rect.Width()*0.05; i <= obj_rect.Width(); i += obj_rect.Width() * 0.1) {
		CString str;
		str.Format(_T("%.3f"), x);
		if (x < 0.01 || x > -0.01) {
			str.Format(_T("%.1f"), x);
		}
		grr->DrawString(str, -1, &font, Gdiplus::PointF(i, 0), &brText);
		x += 0.1 * e_max;
	}

	//рисовка из буфера
	gr.DrawImage(&myBitmap, 0, 0, obj_rect.Width(), obj_rect.Height());
	delete grr;//очистка памяти
}


//СДЕЛАТЬ РИСОВКУ!!!!
void CPotentialpitDlg::DrawWave(std::vector<long double> wave, std::vector<long double> z)
{
	//рисование модели
	CDC* dc = GetDlgItem(IDC_pits)->GetDC();
	HDC hdc = *dc;
	Gdiplus::Graphics gr(hdc); //gr будет рисовать на IDC_Object
	CRect obj_rect; //прямоугольник области рисования
	GetDlgItem(IDC_pits)->GetClientRect(&obj_rect); //получаем размеры прямоугольника
	Gdiplus::Bitmap myBitmap(obj_rect.Width(), obj_rect.Height()); //создаем битовый образ

	Gdiplus::Graphics* grr = Gdiplus::Graphics::FromImage(&myBitmap); //создаем дополнительный объект класса для рисования объектов
	grr->SetSmoothingMode(Gdiplus::SmoothingModeHighSpeed); //устанавливаем сглаживание в режиме наилучшего качества

	grr->Clear(Gdiplus::Color::WhiteSmoke);//очистим фон


	Gdiplus::Matrix mtx(1.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f);//создаем матрицу с осью у, направленной вверх

	long double wave_max = *max_element(wave.begin(), wave.end()),
		wave_min = *min_element(wave.begin(), wave.end()),
		//z_max = *max_element(z.begin(), z.end()),
		//z_min = *min_element(z.begin(), z.end());
	z_max =R,
		z_min = -R;

		wave_max = max(abs(wave_max), abs(wave_min));
	//mtx.Translate(obj_rect.Width() *0.0, obj_rect.Height() * 0.5, Gdiplus::MatrixOrderAppend);//перенос начала координат
		mtx.Translate(obj_rect.Width() * 0.5, obj_rect.Height() * 0.5, Gdiplus::MatrixOrderAppend);//перенос начала координат
	mtx.Scale(obj_rect.Width() / (z_max - z_min), obj_rect.Height() / (2*wave_max)); //изменяем масштаб
	grr->SetTransform(&mtx); //применяем преобразования


	/// ВМЕСТО z[0]/1000 подставить 0.03
	//Gdiplus::Pen pGrid(Gdiplus::Color::Color(50, 50, 50), z[0]/1000); //ручка для основной сетки
	Gdiplus::Pen pGrid(Gdiplus::Color::Color(50, 50, 50), 0.0003); //ручка для основной сетки
	//double CC = C / (2 * wave_max);
	//grr->DrawLine(&pGrid, (Gdiplus::REAL)0, (Gdiplus::REAL)-wave_max, (Gdiplus::REAL)0, (Gdiplus::REAL)C*0.015);//вертикальная
	//grr->DrawLine(&pGrid, (Gdiplus::REAL)z_min, (Gdiplus::REAL)C * 0.015, (Gdiplus::REAL)0, (Gdiplus::REAL)C * 0.015); //ось 0

	/// ВМЕСТО z[0]/1000 подставить 0.03
	Gdiplus::Pen pObj(Gdiplus::Color::DarkRed, z[0] / 800); //ручка для рисования объектов
	/// ВМЕСТО z[0]/1000 подставить 0.03
	Gdiplus::Pen pPits(Gdiplus::Color::Black, z[0] / 900); //ручка для рисования объектов

	//Gdiplus::SolidBrush brObj(Gdiplus::Color::Black); //ручка для рисования объектов

	int wave_size = wave.size();
	double step = 1./1000.;
	for (int i = 0; i <= wave_size*(1-step); i+=step*wave_size) {
		//grr->DrawLine(&pObj, (Gdiplus::REAL)z[i], (Gdiplus::REAL)wave[i], (Gdiplus::REAL)z[i + step * wave_size], (Gdiplus::REAL)wave[i + step * wave_size]);
		Gdiplus::Pen pSubG(Gdiplus::Color::Color(200, 200, 300), 0.0005);
		grr->DrawLine(&pSubG, (Gdiplus::REAL)z[i], (Gdiplus::REAL)wave[i] , (Gdiplus::REAL)z[i + step * wave_size] , (Gdiplus::REAL)wave[i + step * wave_size] );
		long double potential = PotentialFunction(z[i], k_potential/10.,C/10.);
		long double potential_next = PotentialFunction(z[i + step * wave_size], k_potential/10.,C/10.);


		Gdiplus::Pen pSubnG(Gdiplus::Color::Color(10, 10, 0), 0.005);
		grr->DrawLine(&pGrid, (Gdiplus::REAL)z[i], (Gdiplus::REAL)potential - wave_max, (Gdiplus::REAL)z[i+ step * wave_size], (Gdiplus::REAL)potential_next - wave_max);
	}



	//рисовка из буфера
	gr.DrawImage(&myBitmap, 0, 0, obj_rect.Width(), obj_rect.Height());
	delete grr;//очистка памяти

}



void CPotentialpitDlg::OnBnClickedOk()
{

	UpdateData(TRUE);
	/*if (k < 0 || k > k_max) {
		Beep(5000, 500);
		MessageBox(L"НЕ ЛОМАЙ ПРОГРАММУ!!!");
		return;
	}
	BeginWaitCursor();*/


	//Создаём ФФ
	phi = CreatePhaseFunction(R, k_max, energy, k_potential,C);


	//Получаем значение конца интервала, значение ФФ на конце интервала и решение ФФ
	long double e_right = energy[energy.size() - 1];
	long double fi_right = phi[phi.size() - 1];
	long double solFi = SolutionPF(k);

	//С помощью метода половинного деления получаем значение энергии, удолетворяющее решению ФФ, и значение ФФ при такой энергии
	long double fik = 0.0;
	ek = MethodDevision(R, solFi, fi_right, e_right, fik, k_potential,C);


	wave = CreateWaveFunction(R, ek, z, fik, k_potential,C);
	//волновая функция
	
	
	UpdateData(FALSE);
	EndWaitCursor();
	DrawPhi(phi, energy); //фи от энергии
	DrawWave(wave, z);
}
