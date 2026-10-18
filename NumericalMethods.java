
import net.objecthunter.exp4j.Expression;
import net.objecthunter.exp4j.ExpressionBuilder;
import java.util.Scanner;

public class NumericalMethods {

    public static void bisection(Expression e, Scanner input) {
        System.out.println("==========Bisection==========");
        System.out.print("Enter a (lower bound): ");
        double a = Double.parseDouble(input.nextLine());
        System.out.print("Enter b (upper bound): ");
        double b = Double.parseDouble(input.nextLine());
        System.out.print("Enter tolerance (press Enter for default 0.0001): ");
        String tolInput = input.nextLine();
        double tolerance = tolInput.isEmpty() ? 0.0001 : Double.parseDouble(tolInput);
        System.out.print("Enter max iterations (press Enter for default 100): ");
        String iterInput = input.nextLine();
        int maxIter = iterInput.isEmpty() ? 100 : Integer.parseInt(iterInput);

        e.setVariable("x", a);
        double fa = e.evaluate();
        e.setVariable("x", b);
        double fb = e.evaluate();

        if (fa * fb > 0) {
            System.out.println("Error: f(a) and f(b) must have opposite signs.");
            return;
        }

        System.out.printf("%-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s%n",
                "Iteration", "a", "b", "c", "f(a)", "f(b)", "f(c)", "Sign", "Update");
        System.out.println("-".repeat(107));

        int iter = 1;
        while (Math.abs(b - a) > tolerance && iter <= maxIter) {

            e.setVariable("x", a);
            fa = e.evaluate();
            e.setVariable("x", b);
            fb = e.evaluate();

            double c = (a + b) / 2;
            e.setVariable("x", c);
            double fc = e.evaluate();

            String sign = (fa * fc < 0) ? "negative" : "positive";
            String update;
            if (fa * fc < 0) {
                update = "b = c";
            } else {
                update = "a = c";
            }

            System.out.printf("%-10d | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10s | %-10s%n",
                    iter, a, b, c, fa, fb, fc, sign, update);

            if (fa * fc < 0) {
                b = c;
            } else {
                a = c;
            }

            iter++;
        }
        System.out.println("\nRoot ≈ " + (a + b) / 2);
    }

    public static void newtonRaphson(Expression e, Scanner input) {
        System.out.println("==========Newton-Raphson==========");
        System.out.print("Enter x0 (initial guess): ");
        double x0 = Double.parseDouble(input.nextLine());
        System.out.print("Enter tolerance (press Enter for default 0.0001): ");
        String tolInput = input.nextLine();
        double tolerance = tolInput.isEmpty() ? 0.0001 : Double.parseDouble(tolInput);
        System.out.print("Enter max iterations (press Enter for default 100): ");
        String iterInput = input.nextLine();
        int maxIter = iterInput.isEmpty() ? 100 : Integer.parseInt(iterInput);

        double h = 0.0001;

        System.out.printf("%-10s | %-10s | %-10s | %-10s | %-10s | %-10s%n",
                "Iteration", "x0", "f(x0)", "f'(x0)", "x1", "Error");
        System.out.println("-".repeat(75));

        int iter = 1;
        while (iter <= maxIter) {
            e.setVariable("x", x0);
            double fx0 = e.evaluate();

            e.setVariable("x", x0 + h);
            double fxh = e.evaluate();
            double fPrime = (fxh - fx0) / h;

            if (Math.abs(fPrime) < 1e-10) {
                System.out.println("Error: Derivative is zero. Method failed.");
                return;
            }

            double x1 = x0 - fx0 / fPrime;
            double error = Math.abs(x1 - x0);

            // first iteration has no error
            if (iter == 1) {
                System.out.printf("%-10d | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10s%n",
                        iter, x0, fx0, fPrime, x1, "-");
            } else {
                System.out.printf("%-10d | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10.6f%n",
                        iter, x0, fx0, fPrime, x1, error);
            }

            if (error < tolerance) {
                x0 = x1;
                break;
            }

            x0 = x1;
            iter++;
        }
        System.out.println("\nRoot ≈ " + x0);
    }

    public static void secant(Expression e, Scanner input) {
        System.out.println("==========Secant==========");
        System.out.print("Enter x0 (first guess): ");
        double x0 = Double.parseDouble(input.nextLine());
        System.out.print("Enter x1 (second guess): ");
        double x1 = Double.parseDouble(input.nextLine());
        System.out.print("Enter tolerance (press Enter for default 0.0001): ");
        String tolInput = input.nextLine();
        double tolerance = tolInput.isEmpty() ? 0.0001 : Double.parseDouble(tolInput);
        System.out.print("Enter max iterations (press Enter for default 100): ");
        String iterInput = input.nextLine();
        int maxIter = iterInput.isEmpty() ? 100 : Integer.parseInt(iterInput);

        if (x0 == x1) {
            System.out.println("Error: x0 and x1 must not be equal.");
            return;
        }

        System.out.printf("%-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s%n",
                "Iteration", "Xn-1", "Xn", "f(Xn-1)", "f(Xn)", "Xn+1", "Error");
        System.out.println("-".repeat(87));

        int iter = 1;
        while (iter <= maxIter) {
            e.setVariable("x", x0);
            double fx0 = e.evaluate();

            e.setVariable("x", x1);
            double fx1 = e.evaluate();

            if (Math.abs(fx1 - fx0) < 1e-10) {
                System.out.println("Error: Division by zero. Method failed.");
                return;
            }

            double x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
            double error = Math.abs(x2 - x1);

            System.out.printf("%-10d | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10.6f%n",
                    iter, x0, x1, fx0, fx1, x2, error);

            if (error < tolerance) {
                x1 = x2;
                break;
            }

            x0 = x1;
            x1 = x2;
            iter++;
        }
        System.out.println("\nRoot ≈ " + x1);
    }

    public static void forwardPoint(Expression e, Scanner input) {
        System.out.println("==========Forward Difference==========");
        System.out.print("Enter x: ");
        double x = Double.parseDouble(input.nextLine());
        System.out.print("Enter h: ");
        double h = Double.parseDouble(input.nextLine());

        if (h == 0) {
            System.out.println("Error: h cannot be zero.");
            return;
        }

        e.setVariable("x", x);
        double fx = e.evaluate();

        e.setVariable("x", x + h);
        double fxh = e.evaluate();

        double result = (fxh - fx) / h;

        System.out.printf("%-15s | %-15s | %-15s | %-15s%n",
                "x", "f(x)", "f(x+h)", "f'(x)");
        System.out.println("-".repeat(67));
        System.out.printf("%-15.6f | %-15.6f | %-15.6f | %-15.6f%n",
                x, fx, fxh, result);

        System.out.println("\nDerivative ≈ " + result);
    }

    public static void points3(Expression e, Scanner input) {
        System.out.println("==========Three-Point Formula==========");
        System.out.print("Enter x: ");
        double x = Double.parseDouble(input.nextLine());
        System.out.print("Enter h: ");
        double h = Double.parseDouble(input.nextLine());

        if (h == 0) {
            System.out.println("Error: h cannot be zero.");
            return;
        }

        e.setVariable("x", x - h);
        double fxmh = e.evaluate();

        e.setVariable("x", x);
        double fx = e.evaluate();

        e.setVariable("x", x + h);
        double fxph = e.evaluate();

        double result = (fxph - fxmh) / (2 * h);

        System.out.printf("%-15s | %-15s | %-15s | %-15s | %-15s%n",
                "x", "f(x-h)", "f(x)", "f(x+h)", "f'(x)");
        System.out.println("-".repeat(83));
        System.out.printf("%-15.6f | %-15.6f | %-15.6f | %-15.6f | %-15.6f%n",
                x, fxmh, fx, fxph, result);

        System.out.println("\nDerivative ≈ " + result);
    }

    public static void points5(Expression e, Scanner input) {
        System.out.println("==========Five-Point Formula==========");
        System.out.print("Enter x: ");
        double x = Double.parseDouble(input.nextLine());
        System.out.print("Enter h: ");
        double h = Double.parseDouble(input.nextLine());

        if (h == 0) {
            System.out.println("Error: h cannot be zero.");
            return;
        }

        e.setVariable("x", x - 2 * h);
        double fx2mh = e.evaluate();

        e.setVariable("x", x - h);
        double fxmh = e.evaluate();

        e.setVariable("x", x);
        double fx = e.evaluate();

        e.setVariable("x", x + h);
        double fxph = e.evaluate();

        e.setVariable("x", x + 2 * h);
        double fx2ph = e.evaluate();

        double result = (-fx2ph + 8 * fxph - 8 * fxmh + fx2mh) / (12 * h);

        System.out.printf("%-15s | %-15s | %-15s | %-15s | %-15s | %-15s | %-15s%n",
                "x", "f(x-2h)", "f(x-h)", "f(x)", "f(x+h)", "f(x+2h)", "f'(x)");
        System.out.println("-".repeat(115));
        System.out.printf("%-15.6f | %-15.6f | %-15.6f | %-15.6f | %-15.6f | %-15.6f | %-15.6f%n",
                x, fx2mh, fxmh, fx, fxph, fx2ph, result);

        System.out.println("\nDerivative ≈ " + result);
    }

    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);
        System.out.println("Supported functions:");
        System.out.println("  log(x)  --> natural log (ln)");
        System.out.println("  log2(x) --> log base 2");
        System.out.println("  sin(x), cos(x), tan(x)");
        System.out.println("  sqrt(x), abs(x)");
        System.out.println("  Use ^ for powers e.g. x^2");
        System.out.println("Enter your formula:");
        String formula = input.nextLine();
        Expression e = new ExpressionBuilder(formula)
                .variable("x")
                .build();
        System.out.println("\nChoose a method:");
        System.out.println("1. Bisection Method");
        System.out.println("2. Newton Raphson Method");
        System.out.println("3. Secant Method");
        System.out.println("4. Forward Point Formula");
        System.out.println("5. Three-Point Formula");
        System.out.println("6. Five-Point Formula");
        System.out.print("Enter your choice: ");
        int choice = Integer.parseInt(input.nextLine());
        if (choice == 1) {
            bisection(e, input);
        } else if (choice == 2) {
            newtonRaphson(e, input);
        } else if (choice == 3) {
            secant(e, input);
        } else if (choice == 4) {
            forwardPoint(e, input);
        } else if (choice == 5) {
            points3(e, input);
        } else if (choice == 6) {
            points5(e, input);
        } else {
            System.out.println("Invalid choice.");
        }
    }
}
