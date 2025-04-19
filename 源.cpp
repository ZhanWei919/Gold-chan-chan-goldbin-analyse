#include <iostream>
#include <vector>
#include <cmath>    // For pow, fabs
#include <limits>   // For numeric_limits
#include <iomanip>  // For setprecision



// --- 数值计算参数 ---
const double TOLERANCE = 1e-9;      // 二分法精度容忍度
const int MAX_ITERATIONS = 1000;    // 二分法最大迭代次数，防止死循环

// --- 辅助函数：计算给定p时选项2的期望收益 ---
// E2(p) = (y / 8) * (1 + p + p^2 + ... + p^7)
double calculate_expected_option2(double y, double p) {
    if (y == 0.0) return 0.0; // 如果y是0，收益总是0
    // 直接计算级数和，避免 p=1 时 (1-p^8)/(1-p) 的除零问题
    double sum_powers = 0.0;
    for (int j = 0; j <= 7; ++j) {
        sum_powers += std::pow(p, j);
    }
    return (y / 8.0) * sum_powers;
}

// --- 函数：使用二分法寻找纳什均衡 p* ---
// 求解方程: calculate_expected_option2(y, p) - x = 0
double find_nash_equilibrium_p(double x, double y) {
    // 定义目标函数 f(p) = E2(p) - x
    auto f = [&](double p) {
        return calculate_expected_option2(y, p) - x;
    };

    double low = 0.0;
    double high = 1.0;

    // 检查边界值，理论上在 (y/8, y) 区间内 f(0) < 0, f(1) > 0
    if (f(low) >= 0) {
        // 如果 f(0) >= 0，意味着 x <= y/8，均衡是 p*=0
        // 这种情况理论上应该在调用此函数前被 main 函数处理掉
        // 但为健壮性保留检查
        return 0.0;
    }
    if (f(high) <= 0) {
        // 如果 f(1) <= 0，意味着 x >= y，均衡是 p*=1
        // 这种情况理论上应该在调用此函数前被 main 函数处理掉
        return 1.0;
    }

    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double mid = low + (high - low) / 2.0; // 使用更稳健的中点计算
        double f_mid = f(mid);

        if (std::fabs(f_mid) < TOLERANCE || (high - low) / 2.0 < TOLERANCE) {
            // 如果函数值接近0或者区间足够小，找到解
            return mid;
        }

        if (f_mid < 0) {
            // 如果 f(mid) < 0，说明 E2(mid) < x，需要增大 p 才能让 E2 接近 x
            low = mid;
        }
        else {
            // 如果 f(mid) > 0，说明 E2(mid) > x，需要减小 p
            high = mid;
        }
    }

    // 如果达到最大迭代次数仍未收敛（理论上不应发生在此问题中，除非容差过小）
    std::cerr << "Warning: Bisection method did not converge within max iterations." << std::endl;
    return low + (high - low) / 2.0; // 返回当前最佳近似值
}

int main() {
    double x, y;
    double P_STATISTICAL;
    
    std::cout << "请输入统计到的p: ";
    std::cin >> P_STATISTICAL;
    std::cout << "请输入固定收益 X: ";
    std::cin >> x;
    std::cout << "请输入共享总收益 Y: ";
    std::cin >> y;

    // 输入验证 (基础)
    if (x < 0 || y < 0) {
        std::cerr << "错误: 收益 X 和 Y 不能为负数。" << std::endl;
        return 1;
    }

    double p_star; // 纳什均衡概率

    // 1. 处理纯策略纳什均衡的边界情况
    if (y <= 0 || x >= y) { // 如果 y=0 或 x>=y
        p_star = 1.0;
        std::cout << "分析: X >= Y 或 Y <= 0. 纯策略纳什均衡 p* = 1.0 (所有人都应选择选项1)" << std::endl;
    }
    else if (x <= y / 8.0) {
        p_star = 0.0;
        std::cout << "分析: X <= Y/8. 纯策略纳什均衡 p* = 0.0 (所有人都应选择选项2)" << std::endl;
    }
    else {
        // 2. 计算混合策略纳什均衡 p*
        std::cout << "分析: Y/8 < X < Y. 寻找混合策略纳什均衡 p*..." << std::endl;
        p_star = find_nash_equilibrium_p(x, y);
        std::cout << "计算得到的纳什均衡概率 p* ≈ " << std::fixed << std::setprecision(8) << p_star << std::endl;
    }

    std::cout << "您事先统计得到的群体概率 p = " << P_STATISTICAL << std::endl;

    // 3. 比较 p_star 和 P_STATISTICAL 并给出建议
    std::cout << "\n--- 决策建议 ---" << std::endl;
    // 使用一个小容差进行浮点数比较
    if (std::fabs(P_STATISTICAL - p_star) < TOLERANCE) {
        std::cout << "统计概率 p ≈ 纳什均衡 p*. 选择选项1或选项2的期望收益几乎相同。" << std::endl;
        std::cout << "您可以根据个人风险偏好或其他因素决策。" << std::endl;
    }
    else if (P_STATISTICAL > p_star) {
        // 您认为其他人更倾向于选1（比均衡点更倾向）
        // 这使得选项2的期望收益 E2(p) > E2(p*) = X
        std::cout << "因为 统计概率 p > 纳什均衡 p*, 选择选项2 (平分金币) 的期望收益更高。" << std::endl;
        std::cout << "推荐选择: 选项2 (平分 Y 个金币)" << std::endl;
    }
    else { // P_STATISTICAL < p_star
        // 您认为其他人不那么倾向于选1（比均衡点更不倾向）
        // 这使得选项2的期望收益 E2(p) < E2(p*) = X
        std::cout << "因为 统计概率 p < 纳什均衡 p*, 选择选项1 (独享金币) 的期望收益更高。" << std::endl;
        std::cout << "推荐选择: 选项1 (获取 X 个金币)" << std::endl;
    }

    return 0;
}