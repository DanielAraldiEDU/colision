#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <random>

// =========================================================
// 1. ESTRUTURAS MATEMÁTICAS E CONSTANTES
// =========================================================

// Custos definidos no artigo
const float CUSTO_ESFERA = 10.0f;
const float CUSTO_CAIXA = 200.0f;
const float CUSTO_TRIANGULO = 100.0f; // Assumido constante conforme texto

// Contadores Globais para Métricas
long long g_TotalVolumeTests = 0;
long long g_TotalTriangleTests = 0;

struct Vec3
{
	float x, y, z;
	Vec3(float _x = 0, float _y = 0, float _z = 0) : x(_x), y(_y), z(_z) {}
	Vec3 operator+(const Vec3 &v) const { return {x + v.x, y + v.y, z + v.z}; }
	Vec3 operator-(const Vec3 &v) const { return {x - v.x, y - v.y, z - v.z}; }
	Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
	float dot(const Vec3 &v) const { return x * v.x + y * v.y + z * v.z; }
	Vec3 cross(const Vec3 &v) const { return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x}; }
	float length() const { return std::sqrt(x * x + y * y + z * z); }

	// Rotação simples em torno do eixo Z (para posicionamento no círculo)
	Vec3 rotateZ(float angle) const
	{
		float c = std::cos(angle);
		float s = std::sin(angle);
		return {x * c - y * s, x * s + y * c, z};
	}
};

// Matriz de Transformação 4x4 (Simplificada para Rotação/Translação)
struct Matrix4
{
	float m[4][4];
	Matrix4()
	{
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				m[i][j] = (i == j) ? 1.0f : 0.0f;
	}

	static Matrix4 Translate(const Vec3 &v)
	{
		Matrix4 mat;
		mat.m[0][3] = v.x;
		mat.m[1][3] = v.y;
		mat.m[2][3] = v.z;
		return mat;
	}

	Vec3 transformPoint(const Vec3 &p) const
	{
		return Vec3(
				m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3],
				m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3],
				m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3]);
	}
};

struct Triangulo
{
	Vec3 v[3];
	Vec3 centro;
	Triangulo(Vec3 a, Vec3 b, Vec3 c)
	{
		v[0] = a;
		v[1] = b;
		v[2] = c;
		centro = (a + b + c) * (1.0f / 3.0f);
	}
};

// =========================================================
// 2. GERADORES DE GEOMETRIA (Torus, Pedra, Cubo)
// =========================================================

std::vector<Triangulo> GerarCubo()
{
	// Cubo simples ~12 a 14 triangulos
	std::vector<Triangulo> tris;
	Vec3 p[8] = {
			{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}};
	int indices[] = {
			0, 1, 2, 0, 2, 3, // Back
			4, 5, 6, 4, 6, 7, // Front
			0, 4, 7, 0, 7, 3, // Left
			1, 5, 6, 1, 6, 2, // Right
			3, 2, 6, 3, 6, 7, // Top
			0, 1, 5, 0, 5, 4	// Bottom
	};
	for (int i = 0; i < 36; i += 3)
	{
		tris.emplace_back(p[indices[i]], p[indices[i + 1]], p[indices[i + 2]]);
	}
	return tris; // 12 triângulos (próximo de 14 citado)
}

std::vector<Triangulo> GerarTorus(int rings = 30, int sides = 30, float r1 = 1.0f, float r2 = 3.0f)
{
	// Torus ~1800 triângulos (30 * 30 * 2 = 1800)
	std::vector<Triangulo> tris;
	std::vector<Vec3> verts;

	for (int i = 0; i <= rings; ++i)
	{
		float theta = i * 2.0f * 3.14159f / rings;
		for (int j = 0; j <= sides; ++j)
		{
			float phi = j * 2.0f * 3.14159f / sides;
			float x = (r2 + r1 * std::cos(phi)) * std::cos(theta);
			float y = (r2 + r1 * std::cos(phi)) * std::sin(theta);
			float z = r1 * std::sin(phi);
			verts.push_back({x, y, z});
		}
	}

	for (int i = 0; i < rings; ++i)
	{
		for (int j = 0; j < sides; ++j)
		{
			int current = i * (sides + 1) + j;
			int next = current + sides + 1;
			tris.emplace_back(verts[current], verts[next], verts[current + 1]);
			tris.emplace_back(verts[current + 1], verts[next], verts[next + 1]);
		}
	}
	return tris;
}

std::vector<Triangulo> GerarPedra(int stacks = 20, int slices = 19)
{
	// Pedra (Esfera deformada) ~760 triângulos (20 * 19 * 2 = 760)
	std::vector<Triangulo> tris;
	std::vector<Vec3> verts;

	// Gerar esfera
	for (int i = 0; i <= stacks; ++i)
	{
		float lat = 3.14159f * (-0.5f + (float)i / stacks);
		float z = std::sin(lat);
		float zr = std::cos(lat);

		for (int j = 0; j <= slices; ++j)
		{
			float lng = 2 * 3.14159f * (float)j / slices;
			float x = zr * std::cos(lng);
			float y = zr * std::sin(lng);

			// Adiciona ruído para parecer "Pedra"
			float noise = 1.0f + ((rand() % 100) / 500.0f);
			verts.push_back(Vec3(x * noise, y * noise, z * noise));
		}
	}

	for (int i = 0; i < stacks; ++i)
	{
		for (int j = 0; j < slices; ++j)
		{
			int current = i * (slices + 1) + j;
			int next = current + slices + 1;
			tris.emplace_back(verts[current], verts[current + 1], verts[next]);
			tris.emplace_back(verts[next], verts[current + 1], verts[next + 1]);
		}
	}
	return tris;
}

// =========================================================
// 3. HIERARQUIA DE VOLUMES (BVH)
// =========================================================

enum TipoVolume
{
	ESFERA,
	AABB,
	OBB
};

struct Volume
{
	TipoVolume tipo;
	Vec3 centro;	 // Centro local
	Vec3 extensao; // AABB: half-extents, Esfera: x=raio

	Volume(TipoVolume t) : tipo(t), centro(0, 0, 0), extensao(0, 0, 0) {}
};

struct Nodo
{
	std::unique_ptr<Volume> volume;
	std::vector<Triangulo> triangulos;
	std::unique_ptr<Nodo> esq;
	std::unique_ptr<Nodo> dir;
	bool folha;

	Nodo() : folha(false) {}
};

// Funções de Construção (Simplificadas para brevidade, lógica idêntica à anterior)
std::unique_ptr<Volume> CriarVolume(const std::vector<Triangulo> &tris, TipoVolume tipo)
{
	auto vol = std::make_unique<Volume>(tipo);
	if (tris.empty())
		return vol;

	Vec3 minV(1e9, 1e9, 1e9), maxV(-1e9, -1e9, -1e9);
	Vec3 soma(0, 0, 0);

	for (auto &t : tris)
	{
		for (auto &v : t.v)
		{
			soma = soma + v;
			if (v.x < minV.x)
				minV.x = v.x;
			if (v.x > maxV.x)
				maxV.x = v.x;
			if (v.y < minV.y)
				minV.y = v.y;
			if (v.y > maxV.y)
				maxV.y = v.y;
			if (v.z < minV.z)
				minV.z = v.z;
			if (v.z > maxV.z)
				maxV.z = v.z;
		}
	}
	int numV = tris.size() * 3;
	vol->centro = soma * (1.0f / numV);

	if (tipo == ESFERA)
	{
		float maxDistSq = 0;
		for (auto &t : tris)
			for (auto &v : t.v)
			{
				float d = (v - vol->centro).dot(v - vol->centro);
				if (d > maxDistSq)
					maxDistSq = d;
			}
		vol->extensao.x = std::sqrt(maxDistSq); // Raio
	}
	else
	{
		// AABB e OBB (Simplificado como AABB local para OBB neste teste)
		vol->centro = (minV + maxV) * 0.5f;
		vol->extensao = (maxV - minV) * 0.5f;
	}
	return vol;
}

std::unique_ptr<Nodo> ConstruirBVH(std::vector<Triangulo> &tris, TipoVolume tipo, int depth = 0)
{
	auto nodo = std::make_unique<Nodo>();
	nodo->volume = CriarVolume(tris, tipo);

	// Critério de parada: 4 níveis para cubo, 10 para pedra, 11 para torus
	// Ou limiar de triângulos. Usaremos limiar pequeno para forçar profundidade.
	if (tris.size() <= 2 || depth > 12)
	{
		nodo->folha = true;
		nodo->triangulos = tris;
		return nodo;
	}

	// Subdivisão Espacial Recursiva
	int eixo = depth % 3;
	std::sort(tris.begin(), tris.end(), [eixo](const Triangulo &a, const Triangulo &b)
						{ return (eixo == 0) ? a.centro.x < b.centro.x : (eixo == 1) ? a.centro.y < b.centro.y
																																				 : a.centro.z < b.centro.z; });

	size_t mid = tris.size() / 2;
	std::vector<Triangulo> tEsq(tris.begin(), tris.begin() + mid);
	std::vector<Triangulo> tDir(tris.begin() + mid, tris.end());

	nodo->esq = ConstruirBVH(tEsq, tipo, depth + 1);
	nodo->dir = ConstruirBVH(tDir, tipo, depth + 1);
	return nodo;
}

class ArvoreBVH
{
public:
	std::unique_ptr<Nodo> raiz;
	TipoVolume tipo;
	ArvoreBVH(std::vector<Triangulo> tris, TipoVolume t) : tipo(t)
	{
		raiz = ConstruirBVH(tris, t);
	}
};

// =========================================================
// 4. LÓGICA DE DETECÇÃO (Com Transformações)
// =========================================================

bool TesteTriangulos(const std::vector<Triangulo> &tA, const std::vector<Triangulo> &tB)
{
	g_TotalTriangleTests++; // Métrica T
	// Simulação de custo: Verifica bounding box dos triângulos para evitar N*M
	return true; // Assume colisão se chegou na folha (Simplificação para métrica)
}

// Verifica sobreposição considerando que os volumes estão em coordenadas locais
// e devem ser transformados para o mundo pelas matrizes matA e matB.
bool TesteVolume(const Volume *vA, const Matrix4 &matA, const Volume *vB, const Matrix4 &matB)
{
	g_TotalVolumeTests++; // Métrica V

	// Transforma centros para o World Space
	Vec3 cA = matA.transformPoint(vA->centro);
	Vec3 cB = matB.transformPoint(vB->centro);

	if (vA->tipo == ESFERA)
	{
		float rA = vA->extensao.x;
		float rB = vB->extensao.x;
		float distSq = (cA - cB).dot(cA - cB);
		return distSq <= (rA + rB) * (rA + rB);
	}
	else
	{
		// AABB e OBB (Tratados como OBB pois têm rotação)
		// SAT Simplificado (Apenas teste de distância radial para performance da demo)
		// Numa implementação real de OBB, testaríamos os 15 eixos.
		// Aqui usamos esferas envolventes das caixas para garantir que o contador V seja incrementado
		// mas mantendo a lógica simples.
		float rA = vA->extensao.length();
		float rB = vB->extensao.length();
		return (cA - cB).length() <= (rA + rB);
	}
}

bool TesteRecursivo(Nodo *nA, const Matrix4 &mA, Nodo *nB, const Matrix4 &mB)
{
	// 1. Teste de Volumes
	if (!TesteVolume(nA->volume.get(), mA, nB->volume.get(), mB))
		return false;

	// 2. Descida na Árvore
	if (!nA->folha)
	{
		if (TesteRecursivo(nA->esq.get(), mA, nB, mB))
			return true;
		if (TesteRecursivo(nA->dir.get(), mA, nB, mB))
			return true;
	}
	else if (!nB->folha)
	{
		if (TesteRecursivo(nA, mA, nB->esq.get(), mB))
			return true;
		if (TesteRecursivo(nA, mA, nB->dir.get(), mB))
			return true;
	}
	else
	{
		// 3. Teste de Triângulos
		return TesteTriangulos(nA->triangulos, nB->triangulos);
	}
	return false;
}

// =========================================================
// 5. SIMULAÇÃO E CENÁRIO DE TESTE
// =========================================================

struct CorpoRigido
{
	Vec3 posicao;
	Vec3 velocidade;
	Matrix4 transform;
	ArvoreBVH *arvore; // Referência para geometria compartilhada

	void Update(float dt)
	{
		posicao = posicao + velocidade * dt;
		transform = Matrix4::Translate(posicao); // Atualiza matriz
	}
};

void ExecutarTeste(int nObjetos, TipoVolume tipo, std::string nomeTipo)
{
	// Reset contadores
	g_TotalVolumeTests = 0;
	g_TotalTriangleTests = 0;

	// 1. Preparar Geometria
	std::vector<Triangulo> geoTorus = GerarTorus();
	std::vector<Triangulo> geoPedra = GerarPedra();
	std::vector<Triangulo> geoCubo = GerarCubo();

	ArvoreBVH bvhTorus(geoTorus, tipo);
	ArvoreBVH bvhPedra(geoPedra, tipo);
	ArvoreBVH bvhCubo(geoCubo, tipo);

	// 2. Instanciar Corpos
	std::vector<CorpoRigido> corpos;
	float raioCirculo = 30.0f;

	for (int i = 0; i < nObjetos; i++)
	{
		CorpoRigido corpo;
		// Distribui tipos de objetos alternadamente
		if (i % 3 == 0)
			corpo.arvore = &bvhTorus;
		else if (i % 3 == 1)
			corpo.arvore = &bvhPedra;
		else
			corpo.arvore = &bvhCubo;

		// Posiciona em círculo
		float angulo = (2.0f * 3.14159f * i) / nObjetos;
		corpo.posicao = Vec3(std::cos(angulo) * raioCirculo, std::sin(angulo) * raioCirculo, 0);

		// Velocidade em direção ao centro
		corpo.velocidade = corpo.posicao * -0.05f;
		corpo.Update(0); // Setup inicial matrix
		corpos.push_back(corpo);
	}

	// 3. Loop de Simulação (600 quadros)
	for (int quadro = 0; quadro < 600; quadro++)
	{
		// Atualiza física
		for (auto &c : corpos)
			c.Update(0.1f);

		// Detecção de Colisão Broadphase (N^2)
		for (int i = 0; i < nObjetos; i++)
		{
			for (int j = i + 1; j < nObjetos; j++)
			{
				TesteRecursivo(corpos[i].arvore->raiz.get(), corpos[i].transform,
											 corpos[j].arvore->raiz.get(), corpos[j].transform);
			}
		}
	}

	// 4. Resultados e Equação de Custo
	// Ec = V * Cv + T * Ct
	float C_V = (tipo == ESFERA) ? CUSTO_ESFERA : CUSTO_CAIXA;
	double Ec = (double)g_TotalVolumeTests * C_V + (double)g_TotalTriangleTests * CUSTO_TRIANGULO;

	std::cout << std::left << std::setw(15) << nomeTipo
						<< "| N=" << std::setw(3) << nObjetos
						<< "| Testes Volume (V): " << std::setw(10) << g_TotalVolumeTests
						<< "| Testes Triângulo (T): " << std::setw(8) << g_TotalTriangleTests
						<< "| Custo Total (Ec): " << (long long)Ec << std::endl;
}

int main()
{
	std::srand(1234); // Seed fixa para reprodutibilidade

	std::cout << "=== ANALISE DE DESEMPENHO (Replica Fig 2) ===\n";
	std::cout << "Tipos de Volume: Esfera, AABB, OBB\n";
	std::cout << "Quantidades: 2, 4, 8, 16, 32, 64 objetos\n";
	std::cout << "Objetos: Torus(1800), Pedra(760), Cubo(12-14)\n\n";

	// Executa bateria de testes conforme Seção 4
	int qtds[] = {2, 4, 8, 16, 32, 64};

	for (int n : qtds)
	{
		std::cout << "--- CENA COM " << n << " OBJETOS ---\n";
		ExecutarTeste(n, ESFERA, "Esfera");
		ExecutarTeste(n, AABB, "AABB");
		ExecutarTeste(n, OBB, "OBB"); // Nota: OBB custo computacional é maior na prática
		std::cout << "\n";
	}

	return 0;
}