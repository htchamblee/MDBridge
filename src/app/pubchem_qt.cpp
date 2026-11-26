#include "pubchem_qt.hpp"

#include <cmath>
#include <memory>
#include <string>

#include <QByteArray>
#include <QCoreApplication>
#include <QEventLoop>
#include <QNetworkAccessManager>
#include <QNetworkProxy>
#include <QNetworkProxyFactory>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QString>
#include <QTimer>
#include <QUrl>
#include <QVariant>

namespace mdbridge::app {

namespace {

HttpReply qtGet(const std::string& url, const HttpHeaders& headers, double timeoutSeconds) {
    std::unique_ptr<QCoreApplication> own;
    if (QCoreApplication::instance() == nullptr) {
        static int argc = 1;
        static char name[] = "mdbridge";
        static char* argv[] = {name, nullptr};
        own = std::make_unique<QCoreApplication>(argc, argv);
    }
    QNetworkProxyFactory::setUseSystemConfiguration(true);
    QNetworkAccessManager manager;
    const QUrl target(QString::fromStdString(url));
    const QString host = target.host();
    if (host == QLatin1String("127.0.0.1") || host == QLatin1String("localhost") || host == QLatin1String("::1")) {
        manager.setProxy(QNetworkProxy(QNetworkProxy::NoProxy));  // the loopback is never proxied
    }
    QNetworkRequest request(target);
    for (const auto& [name, value] : headers) {
        request.setRawHeader(QByteArray::fromStdString(name), QByteArray::fromStdString(value));
    }
    request.setAttribute(QNetworkRequest::RedirectPolicyAttribute, QNetworkRequest::NoLessSafeRedirectPolicy);

    QNetworkReply* reply = manager.get(request);
    QEventLoop loop;
    QTimer timer;
    timer.setSingleShot(true);
    bool timedOut = false;
    QObject::connect(&timer, &QTimer::timeout, &loop, [&] {
        timedOut = true;
        reply->abort();
    });
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    const double ms = std::ceil(timeoutSeconds * 1000.0);
    timer.start(ms >= 2147483647.0 ? 2147483647 : ms < 1.0 ? 1 : static_cast<int>(ms));
    if (!reply->isFinished()) {
        loop.exec();
    }
    timer.stop();

    HttpReply out;
    const QVariant status = reply->attribute(QNetworkRequest::HttpStatusCodeAttribute);
    if (timedOut) {
        out.outcome = HttpReply::Outcome::TimedOut;
    } else if (status.isValid()) {
        out.outcome = HttpReply::Outcome::Answered;
        out.status = status.toInt();
        out.reason = reply->attribute(QNetworkRequest::HttpReasonPhraseAttribute).toString().toStdString();
        out.body = reply->readAll().toStdString();
    } else {
        out.outcome = reply->error() == QNetworkReply::TimeoutError ? HttpReply::Outcome::TimedOut
                                                                     : HttpReply::Outcome::Failed;
        out.error = reply->errorString().toStdString();
    }
    return out;  // the reply goes with the manager
}

}  // namespace

NetworkBackend qtNetworkBackend() { return NetworkBackend{&qtGet}; }

}  // namespace mdbridge::app
