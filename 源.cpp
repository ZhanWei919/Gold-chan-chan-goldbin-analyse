#include <iostream>
#include <vector>
#include <cmath>    // For pow, fabs
#include <limits>   // For numeric_limits
#include <iomanip>  // For setprecision



// --- ��ֵ������� ---
const double TOLERANCE = 1e-9;      // ���ַ��������̶�
const int MAX_ITERATIONS = 1000;    // ���ַ���������������ֹ��ѭ��

// --- �����������������pʱѡ��2���������� ---
// E2(p) = (y / 8) * (1 + p + p^2 + ... + p^7)
double calculate_expected_option2(double y, double p) {
    if (y == 0.0) return 0.0; // ���y��0����������0
    // ֱ�Ӽ��㼶���ͣ����� p=1 ʱ (1-p^8)/(1-p) �ĳ�������
    double sum_powers = 0.0;
    for (int j = 0; j <= 7; ++j) {
        sum_powers += std::pow(p, j);
    }
    return (y / 8.0) * sum_powers;
}

// --- ������ʹ�ö��ַ�Ѱ����ʲ���� p* ---
// ��ⷽ��: calculate_expected_option2(y, p) - x = 0
double find_nash_equilibrium_p(double x, double y) {
    // ����Ŀ�꺯�� f(p) = E2(p) - x
    auto f = [&](double p) {
        return calculate_expected_option2(y, p) - x;
    };

    double low = 0.0;
    double high = 1.0;

    // ���߽�ֵ���������� (y/8, y) ������ f(0) < 0, f(1) > 0
    if (f(low) >= 0) {
        // ��� f(0) >= 0����ζ�� x <= y/8�������� p*=0
        // �������������Ӧ���ڵ��ô˺���ǰ�� main ���������
        // ��Ϊ��׳�Ա������
        return 0.0;
    }
    if (f(high) <= 0) {
        // ��� f(1) <= 0����ζ�� x >= y�������� p*=1
        // �������������Ӧ���ڵ��ô˺���ǰ�� main ���������
        return 1.0;
    }

    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double mid = low + (high - low) / 2.0; // ʹ�ø��Ƚ����е����
        double f_mid = f(mid);

        if (std::fabs(f_mid) < TOLERANCE || (high - low) / 2.0 < TOLERANCE) {
            // �������ֵ�ӽ�0���������㹻С���ҵ���
            return mid;
        }

        if (f_mid < 0) {
            // ��� f(mid) < 0��˵�� E2(mid) < x����Ҫ���� p ������ E2 �ӽ� x
            low = mid;
        }
        else {
            // ��� f(mid) > 0��˵�� E2(mid) > x����Ҫ��С p
            high = mid;
        }
    }

    // ����ﵽ������������δ�����������ϲ�Ӧ�����ڴ������У������ݲ��С��
    std::cerr << "Warning: Bisection method did not converge within max iterations." << std::endl;
    return low + (high - low) / 2.0; // ���ص�ǰ��ѽ���ֵ
}

int main() {
    double x, y;
    double P_STATISTICAL;
    
    std::cout << "������ͳ�Ƶ���p: ";
    std::cin >> P_STATISTICAL;
    std::cout << "������̶����� X: ";
    std::cin >> x;
    std::cout << "�����빲�������� Y: ";
    std::cin >> y;

    // ������֤ (����)
    if (x < 0 || y < 0) {
        std::cerr << "����: ���� X �� Y ����Ϊ������" << std::endl;
        return 1;
    }

    double p_star; // ��ʲ�������

    // 1. ����������ʲ����ı߽����
    if (y <= 0 || x >= y) { // ��� y=0 �� x>=y
        p_star = 1.0;
        std::cout << "����: X >= Y �� Y <= 0. ��������ʲ���� p* = 1.0 (�����˶�Ӧѡ��ѡ��1)" << std::endl;
    }
    else if (x <= y / 8.0) {
        p_star = 0.0;
        std::cout << "����: X <= Y/8. ��������ʲ���� p* = 0.0 (�����˶�Ӧѡ��ѡ��2)" << std::endl;
    }
    else {
        // 2. �����ϲ�����ʲ���� p*
        std::cout << "����: Y/8 < X < Y. Ѱ�һ�ϲ�����ʲ���� p*..." << std::endl;
        p_star = find_nash_equilibrium_p(x, y);
        std::cout << "����õ�����ʲ������� p* �� " << std::fixed << std::setprecision(8) << p_star << std::endl;
    }

    std::cout << "������ͳ�Ƶõ���Ⱥ����� p = " << P_STATISTICAL << std::endl;

    // 3. �Ƚ� p_star �� P_STATISTICAL ����������
    std::cout << "\n--- ���߽��� ---" << std::endl;
    // ʹ��һ��С�ݲ���и������Ƚ�
    if (std::fabs(P_STATISTICAL - p_star) < TOLERANCE) {
        std::cout << "ͳ�Ƹ��� p �� ��ʲ���� p*. ѡ��ѡ��1��ѡ��2���������漸����ͬ��" << std::endl;
        std::cout << "�����Ը��ݸ��˷���ƫ�û��������ؾ��ߡ�" << std::endl;
    }
    else if (P_STATISTICAL > p_star) {
        // ����Ϊ�����˸�������ѡ1���Ⱦ���������
        // ��ʹ��ѡ��2���������� E2(p) > E2(p*) = X
        std::cout << "��Ϊ ͳ�Ƹ��� p > ��ʲ���� p*, ѡ��ѡ��2 (ƽ�ֽ��) ������������ߡ�" << std::endl;
        std::cout << "�Ƽ�ѡ��: ѡ��2 (ƽ�� Y �����)" << std::endl;
    }
    else { // P_STATISTICAL < p_star
        // ����Ϊ�����˲���ô������ѡ1���Ⱦ�����������
        // ��ʹ��ѡ��2���������� E2(p) < E2(p*) = X
        std::cout << "��Ϊ ͳ�Ƹ��� p < ��ʲ���� p*, ѡ��ѡ��1 (������) ������������ߡ�" << std::endl;
        std::cout << "�Ƽ�ѡ��: ѡ��1 (��ȡ X �����)" << std::endl;
    }

    return 0;
}